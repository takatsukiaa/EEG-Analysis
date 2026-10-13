function results = get_eeg_windowed_mtlcca_artifact_removal(eeg_data, Fs, options)
%GET_EEG_WINDOWED_MTLCCA_ARTIFACT_REMOVAL Conservative windowed MTLCCA.
%   RESULTS = GET_EEG_WINDOWED_MTLCCA_ARTIFACT_REMOVAL(EEG_DATA, FS)
%   separates lagged copies of multichannel EEG jointly with SSQCOR
%   multiset CCA. Components with low between-lag predictability and
%   corroborating spectral/topographic/EMG evidence are softly attenuated.
%
%   The decomposition follows the multi-dataset CCA construction described
%   by Shen et al. (EMBC 2013, doi:10.1109/EMBC.2013.6611116), using the
%   deflationary SSQCOR objective of Li et al. (IEEE TSP 2009,
%   doi:10.1109/TSP.2009.2021636). Soft gains and preservation guardrails
%   are conservative extensions; automatic processing never hard-zeros a
%   component.
%
%   Important options (name=value):
%       WindowDurationSec             default 8
%       OverlapFraction               default 0.75
%       LagSamples                    [] => 10 lags spanning 0-40 ms
%       NumLagDatasets                default 10
%       MaxLagDurationSec             default 0.040
%       EMGReference                  optional synchronized vector
%       PeripheralChannelNames        labels used for topographic evidence
%       MaxCandidateFraction          default 0.25
%       MinComponentGain              default 0.25 (-12.0 dB)
%       ComponentGainsOverride        [] or numeric/cell gains in [0,1]
%       ApplyPreservationGuard        default true
%
%   During a generalized tonic-clonic seizure, true fast ictal activity can
%   resemble EMG. Always compare the output with the untouched recording.

    arguments
        eeg_data (1,1) struct
        Fs (1,1) double {mustBeFinite,mustBePositive}
        options.WindowDurationSec (1,1) double ...
            {mustBeFinite,mustBePositive} = 8
        options.OverlapFraction (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.75
        options.LagSamples = []
        options.NumLagDatasets (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBeGreaterThanOrEqual(options.NumLagDatasets,2)} = 10
        options.MaxLagDurationSec (1,1) double ...
            {mustBeFinite,mustBePositive} = 0.040
        options.CovarianceRegularization (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 1e-6
        options.MCCAMaxSweeps (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 200
        options.MCCATolerance (1,1) double ...
            {mustBeFinite,mustBePositive} = 1e-8
        options.MCCAInitializations (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 3
        options.MaxCandidateFraction (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.25
        options.MinComponentGain (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.25
        options.EMGReference = []
        options.PeripheralChannelNames = ["Fp1","Fp2","F7","F8", ...
            "T3","T4","T5","T6","Fpz","A1","A2"]
        options.ComponentGainsOverride = []
        options.ApplyPreservationGuard (1,1) logical = true
        options.MinLowBandCorrelation (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.98
        options.MaxLowBandPowerChangeDb (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.50
        options.MinSpatialCovarianceCorrelation (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.95
        options.RobustClipMAD (1,1) double ...
            {mustBeFinite,mustBePositive} = 12
    end

    validate_options(options);
    require_signal_toolbox();
    [eeg_matrix, channel_names] = eeg_structure_to_matrix(eeg_data);
    [num_samples, num_channels] = size(eeg_matrix);
    if num_channels < 2
        error('get_eeg_windowed_mtlcca_artifact_removal:NotEnoughEEGChannels', ...
            'At least two EEG channels are required.');
    end

    lag_samples = resolve_lag_samples(options, Fs);
    maximum_lag = max(lag_samples);
    requested_window_samples = max(2,round(options.WindowDurationSec*Fs));
    window_samples = min(num_samples,requested_window_samples);
    if window_samples - maximum_lag <= 2*num_channels
        error('get_eeg_windowed_mtlcca_artifact_removal:SignalTooShort', ...
            ['Each aligned window must contain more than twice the number ' ...
             'of EEG channels after removing the maximum lag.']);
    end

    hop_samples = max(1,round(window_samples*(1-options.OverlapFraction)));
    window_starts = make_window_starts(num_samples,window_samples,hop_samples);
    num_windows = numel(window_starts);
    gain_overrides = resolve_gain_overrides( ...
        options.ComponentGainsOverride,num_windows,num_channels);
    emg_reference = resolve_emg_reference(options.EMGReference,num_samples);
    peripheral_mask = ismember(lower(string(channel_names)), ...
        lower(string(options.PeripheralChannelNames)));

    cleaned_accumulator = zeros(size(eeg_matrix));
    weight_accumulator = zeros(num_samples,1);
    window_ends = zeros(num_windows,1);
    between_lag_scores = cell(num_windows,1);
    component_features = cell(num_windows,1);
    component_gains = cell(num_windows,1);
    component_unmixing = cell(num_windows,1);
    component_mixing = cell(num_windows,1);
    lag_unmixing = cell(num_windows,1);
    effective_eeg_rank = zeros(num_windows,1);
    mcca_objective = cell(num_windows,1);
    mcca_sweeps = cell(num_windows,1);
    reconstruction_error = zeros(num_windows,1);
    guard_strength = zeros(num_windows,1);
    preservation_before_guard = repmat(empty_preservation_metrics(),num_windows,1);
    preservation_after_guard = repmat(empty_preservation_metrics(),num_windows,1);
    skipped_reason = strings(num_windows,1);

    for window_index = 1:num_windows
        first_sample = window_starts(window_index);
        last_sample = min(first_sample+window_samples-1,num_samples);
        indices = first_sample:last_sample;
        window_ends(window_index) = last_sample;
        eeg_window = eeg_matrix(indices,:);
        if isempty(emg_reference)
            emg_window = [];
        else
            emg_window = emg_reference(indices);
        end

        [cleaned_window,diagnostics] = clean_mtlcca_window( ...
            eeg_window,emg_window,Fs,lag_samples,peripheral_mask, ...
            gain_overrides{window_index},options);
        weight = overlap_weight(numel(indices));
        cleaned_accumulator(indices,:) = cleaned_accumulator(indices,:) + ...
            cleaned_window.*weight;
        weight_accumulator(indices) = weight_accumulator(indices)+weight;

        between_lag_scores{window_index} = diagnostics.between_lag_scores;
        component_features{window_index} = diagnostics.component_features;
        component_gains{window_index} = diagnostics.component_gains;
        component_unmixing{window_index} = diagnostics.component_unmixing;
        component_mixing{window_index} = diagnostics.component_mixing;
        lag_unmixing{window_index} = diagnostics.lag_unmixing;
        effective_eeg_rank(window_index) = diagnostics.effective_eeg_rank;
        mcca_objective{window_index} = diagnostics.mcca_objective;
        mcca_sweeps{window_index} = diagnostics.mcca_sweeps;
        reconstruction_error(window_index) = diagnostics.reconstruction_error;
        guard_strength(window_index) = diagnostics.guard_strength;
        preservation_before_guard(window_index) = diagnostics.preservation_before_guard;
        preservation_after_guard(window_index) = diagnostics.preservation_after_guard;
        skipped_reason(window_index) = diagnostics.skipped_reason;
    end

    if any(weight_accumulator<=0)
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidWindowCoverage', ...
            'Overlap-add windows did not cover every output sample.');
    end
    cleaned_matrix = cleaned_accumulator./weight_accumulator;
    artifact_matrix = eeg_matrix-cleaned_matrix;

    results = struct();
    results.cleaned_eeg = matrix_to_eeg_structure(cleaned_matrix,channel_names);
    results.artifact_estimate = matrix_to_eeg_structure(artifact_matrix,channel_names);
    results.channel_names = channel_names;
    results.window_start_samples = window_starts(:);
    results.window_end_samples = window_ends;
    results.window_start_sec = (window_starts(:)-1)/Fs;
    results.window_end_sec = window_ends/Fs;
    results.lag_samples = lag_samples;
    results.lag_sec = lag_samples/Fs;
    results.between_lag_scores = between_lag_scores;
    results.component_features = component_features;
    results.component_gains = component_gains;
    results.component_unmixing = component_unmixing;
    results.component_mixing = component_mixing;
    results.lag_unmixing = lag_unmixing;
    results.effective_eeg_rank = effective_eeg_rank;
    results.mcca_objective = mcca_objective;
    results.mcca_sweeps = mcca_sweeps;
    results.reconstruction_error = reconstruction_error;
    results.guard_strength = guard_strength;
    results.preservation_before_guard = preservation_before_guard;
    results.preservation_after_guard = preservation_after_guard;
    results.skipped_reason = skipped_reason;
    results.metrics = calculate_metrics( ...
        eeg_matrix,cleaned_matrix,artifact_matrix,Fs,emg_reference,channel_names);
    results.settings = struct( ...
        'method','Windowed exact SSQCOR multi-time-lag CCA with conservative soft attenuation', ...
        'window_duration_sec',options.WindowDurationSec, ...
        'window_samples',window_samples, ...
        'overlap_fraction',options.OverlapFraction, ...
        'hop_samples',hop_samples, ...
        'lag_samples',lag_samples, ...
        'lag_sec',lag_samples/Fs, ...
        'covariance_regularization',options.CovarianceRegularization, ...
        'mcca_max_sweeps',options.MCCAMaxSweeps, ...
        'mcca_tolerance',options.MCCATolerance, ...
        'mcca_initializations',options.MCCAInitializations, ...
        'max_candidate_fraction',options.MaxCandidateFraction, ...
        'minimum_component_gain',options.MinComponentGain, ...
        'automatic_selection',isempty(options.ComponentGainsOverride), ...
        'apply_preservation_guard',options.ApplyPreservationGuard, ...
        'minimum_low_band_correlation',options.MinLowBandCorrelation, ...
        'maximum_low_band_power_change_db',options.MaxLowBandPowerChangeDb, ...
        'minimum_spatial_covariance_correlation',options.MinSpatialCovarianceCorrelation, ...
        'robust_clip_mad',options.RobustClipMAD, ...
        'has_emg_reference',~isempty(emg_reference), ...
        'peripheral_channel_names',{channel_names(peripheral_mask)});
    results.notes = [ ...
        "Automatic processing never hard-zeros a component.", ...
        "EMG is supporting evidence only and is not an MCCA input.", ...
        "Fast ictal EEG can resemble muscle; compare with the untouched recording."];
end

function validate_options(options)
    if options.OverlapFraction>=1
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidOverlap', ...
            'OverlapFraction must be in [0,1).');
    end
    bounded = {'CovarianceRegularization','MaxCandidateFraction', ...
        'MinComponentGain','MinLowBandCorrelation', ...
        'MinSpatialCovarianceCorrelation'};
    for k=1:numel(bounded)
        value=options.(bounded{k});
        if value>1
            error('get_eeg_windowed_mtlcca_artifact_removal:InvalidOption', ...
                '%s must be in [0,1].',bounded{k});
        end
    end
end

function require_signal_toolbox()
    required={'butter','filtfilt','pwelch','mscohere'};
    missing=required(cellfun(@(f) exist(f,'file')~=2,required));
    if ~isempty(missing)
        error('get_eeg_windowed_mtlcca_artifact_removal:MissingSignalToolbox', ...
            'Required signal-processing functions are missing: %s.',strjoin(missing,', '));
    end
end

function lags=resolve_lag_samples(options,Fs)
    if isempty(options.LagSamples)
        maximum=max(1,options.MaxLagDurationSec*Fs);
        % Space in continuous time first, then resolve each delay to the
        % nearest sample. At 512 Hz over 40 ms this gives
        % [0 2 5 7 9 11 14 16 18 20].
        lags=unique(round(linspace(0,maximum,options.NumLagDatasets)),'stable');
    else
        lags=double(options.LagSamples(:).');
        if any(~isfinite(lags)) || any(lags<0) || any(lags~=round(lags))
            error('get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples', ...
                'LagSamples must contain finite nonnegative integers.');
        end
        if numel(unique(lags)) ~= numel(lags) || any(diff(lags) <= 0)
            error('get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples', ...
                'LagSamples must be strictly increasing and contain no duplicates.');
        end
        lags=unique(lags,'sorted');
    end
    if numel(lags)<2 || lags(1)~=0
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples', ...
            'LagSamples must contain zero and at least one positive lag.');
    end
end

function emg=resolve_emg_reference(requested,num_samples)
    if isempty(requested)
        emg=[];
        return;
    end
    if ~isnumeric(requested) || ~isreal(requested) || ~isvector(requested)
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidEMGReference', ...
            'EMGReference must be a real numeric vector.');
    end
    emg=double(requested(:));
    if numel(emg)~=num_samples || any(~isfinite(emg))
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidEMGReference', ...
            'EMGReference must be finite and match the EEG sample count.');
    end
    emg=emg-median(emg);
    scale=sqrt(mean(emg.^2));
    if scale<=eps(max(1,max(abs(emg))))
        emg=[];
    else
        emg=emg/scale;
    end
end

function [cleaned,diagnostics]=clean_mtlcca_window( ...
        eeg,emg,Fs,lags,peripheral_mask,gain_override,options)
    diagnostics=empty_window_diagnostics();
    cleaned=eeg;
    robust_eeg=robust_clip_copy(eeg,options.RobustClipMAD);
    maximum_lag=max(lags);
    base=(maximum_lag+1):size(eeg,1);
    num_views=numel(lags);
    raw_views=cell(num_views,1);
    for m=1:num_views
        raw_views{m}=robust_eeg(base-lags(m),:);
    end

    [views,whiteners,dewhiteners,centers,common_rank]= ...
        whiten_views(raw_views,options.CovarianceRegularization);
    diagnostics.effective_eeg_rank=common_rank;
    if common_rank==0
        diagnostics.skipped_reason="Aligned EEG window has zero numerical rank";
        return;
    end

    [lag_weights,objectives,sweeps]=solve_ssqcor_mcca( ...
        views,options.MCCAMaxSweeps,options.MCCATolerance, ...
        options.MCCAInitializations);
    [lag_weights,predictability]=sort_mcca_components(views,lag_weights);
    diagnostics.mcca_objective=objectives;
    diagnostics.mcca_sweeps=sweeps;
    diagnostics.between_lag_scores=predictability;

    zero_view=find(lags==0,1,'first');
    zero_unmixing=whiteners{zero_view}*lag_weights{zero_view};
    zero_mixing=pinv(zero_unmixing);
    center=centers{zero_view};
    centered_full=eeg-center;
    full_scores=centered_full*zero_unmixing;
    rank_projection=full_scores*zero_mixing;
    rank_residual=centered_full-rank_projection;
    no_rejection=full_scores*zero_mixing+rank_residual+center;
    diagnostics.reconstruction_error=norm(no_rejection-eeg,'fro')/max(norm(eeg,'fro'),eps);
    if diagnostics.reconstruction_error>1e-8
        diagnostics.skipped_reason="No-attenuation reconstruction exceeded tolerance";
        return;
    end

    aligned_scores=views{zero_view}*lag_weights{zero_view};
    if isempty(emg)
        aligned_emg=[];
    else
        aligned_emg=emg(base);
    end
    features=extract_component_features(aligned_scores,zero_mixing, ...
        predictability,aligned_emg,Fs,peripheral_mask);
    diagnostics.component_features=features;

    if isempty(gain_override)
        gains=automatic_component_gains(features,options);
    else
        gains=validate_window_gain_override(gain_override,common_rank);
    end
    candidate_clean=reconstruct_with_gains( ...
        full_scores,zero_mixing,rank_residual,center,gains);
    diagnostics.preservation_before_guard=preservation_metrics(eeg,candidate_clean,Fs);

    guard_strength=1;
    if options.ApplyPreservationGuard && ...
            ~passes_preservation(diagnostics.preservation_before_guard,options)
        lower=0; upper=1;
        for iteration=1:16
            trial=(lower+upper)/2;
            trial_gains=1-trial*(1-gains);
            trial_clean=reconstruct_with_gains( ...
                full_scores,zero_mixing,rank_residual,center,trial_gains);
            if passes_preservation(preservation_metrics(eeg,trial_clean,Fs),options)
                lower=trial;
            else
                upper=trial;
            end
        end
        guard_strength=lower;
        gains=1-guard_strength*(1-gains);
    end
    cleaned=reconstruct_with_gains( ...
        full_scores,zero_mixing,rank_residual,center,gains);
    diagnostics.preservation_after_guard=preservation_metrics(eeg,cleaned,Fs);
    diagnostics.guard_strength=guard_strength;
    diagnostics.component_gains=gains;
    diagnostics.component_unmixing=zero_unmixing;
    diagnostics.component_mixing=zero_mixing;
    diagnostics.lag_unmixing=cellfun(@(p,w) p*w, ...
        whiteners,lag_weights,'UniformOutput',false);
    diagnostics.skipped_reason="";

    if any(~isfinite(cleaned),'all')
        error('get_eeg_windowed_mtlcca_artifact_removal:NonfiniteReconstruction', ...
            'MTLCCA reconstruction produced NaN or Inf.');
    end
end

function clipped=robust_clip_copy(matrix,threshold)
    location=median(matrix,1);
    deviation=median(abs(matrix-location),1);
    scale=1.4826*deviation;
    fallback=std(matrix,0,1);
    scale(scale<=eps(max(1,max(abs(matrix),[],1))))= ...
        fallback(scale<=eps(max(1,max(abs(matrix),[],1))));
    scale(scale<=0)=1;
    lower=location-threshold*scale;
    upper=location+threshold*scale;
    clipped=min(max(matrix,lower),upper);
end

function [views,whiteners,dewhiteners,centers,common_rank]= ...
        whiten_views(raw_views,regularization)
    num_views=numel(raw_views);
    eig_vectors=cell(num_views,1);
    eig_values=cell(num_views,1);
    centers=cell(num_views,1);
    ranks=zeros(num_views,1);
    for m=1:num_views
        centers{m}=mean(raw_views{m},1);
        centered=raw_views{m}-centers{m};
        covariance=(centered'*centered)/max(size(centered,1)-1,1);
        covariance=(covariance+covariance')/2;
        [vectors,values]=eig(covariance,'vector');
        [values,order]=sort(real(values),'descend');
        vectors=real(vectors(:,order));
        tolerance=max(size(centered))*eps(max([values(:);1]));
        ranks(m)=sum(values>tolerance);
        eig_vectors{m}=vectors;
        eig_values{m}=values;
    end
    common_rank=min(ranks);
    views=cell(num_views,1);
    whiteners=cell(num_views,1);
    dewhiteners=cell(num_views,1);
    if common_rank==0
        return;
    end
    for m=1:num_views
        vectors=eig_vectors{m}(:,1:common_rank);
        values=eig_values{m}(1:common_rank);
        ridge=regularization*mean(values);
        adjusted=max(values+ridge,eps(max(values)));
        whiteners{m}=vectors*diag(1./sqrt(adjusted));
        dewhiteners{m}=diag(sqrt(adjusted))*vectors';
        centered=raw_views{m}-centers{m};
        views{m}=centered*whiteners{m};
        column_scale=sqrt(sum(views{m}.^2,1)/max(size(views{m},1)-1,1));
        column_scale(column_scale<=eps)=1;
        views{m}=views{m}./column_scale;
        whiteners{m}=whiteners{m}./column_scale;
        dewhiteners{m}=diag(column_scale)*dewhiteners{m};
    end
end

function [weights,objectives,sweep_counts]=solve_ssqcor_mcca( ...
        views,max_sweeps,tolerance,num_initializations)
    num_views=numel(views);
    rank_value=size(views{1},2);
    cross=cell(num_views,num_views);
    sample_denominator=max(size(views{1},1)-1,1);
    for m=1:num_views
        for n=1:num_views
            cross{m,n}=(views{m}'*views{n})/sample_denominator;
        end
    end
    weights=cellfun(@(~) zeros(rank_value,rank_value),views,'UniformOutput',false);
    objectives=zeros(rank_value,1);
    sweep_counts=zeros(rank_value,1);
    for component=1:rank_value
        best_objective=-Inf;
        best_vectors=[];
        best_sweeps=0;
        for initialization=1:num_initializations
            vectors=initial_vectors(cross,weights,component,initialization);
            previous=-Inf;
            for sweep=1:max_sweeps
                for m=1:num_views
                    Q=zeros(rank_value);
                    for n=1:num_views
                        if n==m, continue; end
                        v=vectors{n};
                        Q=Q+cross{m,n}*(v*v')*cross{n,m};
                    end
                    previous_basis=weights{m}(:,1:component-1);
                    projector=eye(rank_value)-previous_basis*previous_basis';
                    Q=projector*((Q+Q')/2)*projector;
                    [candidate,value]=leading_eigenvector(Q);
                    if value<=eps
                        candidate=projected_seed(rank_value,component+m+initialization,previous_basis);
                    end
                    candidate=orthogonalize(candidate,previous_basis);
                    vectors{m}=deterministic_sign(candidate);
                end
                objective=ssqcor_objective(cross,vectors);
                if isfinite(previous) && abs(objective-previous)<= ...
                        tolerance*max(1,abs(previous))
                    break;
                end
                previous=objective;
            end
            if objective>best_objective
                best_objective=objective;
                best_vectors=vectors;
                best_sweeps=sweep;
            end
        end
        for m=1:num_views
            weights{m}(:,component)=best_vectors{m};
        end
        objectives(component)=best_objective;
        sweep_counts(component)=best_sweeps;
    end
end

function vectors=initial_vectors(cross,weights,component,initialization)
    num_views=size(cross,1);
    rank_value=size(cross{1,1},1);
    vectors=cell(num_views,1);
    for m=1:num_views
        previous=weights{m}(:,1:component-1);
        if initialization==1
            seed_index=mod(component+m-2,rank_value)+1;
            seed=zeros(rank_value,1); seed(seed_index)=1;
        elseif initialization==2
            Q=zeros(rank_value);
            for n=1:num_views
                if n~=m, Q=Q+cross{m,n}*cross{n,m}; end
            end
            projector=eye(rank_value)-previous*previous';
            [seed,~]=leading_eigenvector(projector*Q*projector);
        else
            index=(1:rank_value)';
            seed=sin((initialization+component)*index)+ ...
                cos((m+component)*index/2);
        end
        vectors{m}=deterministic_sign(orthogonalize(seed,previous));
    end
end

function vector=projected_seed(dimension,seed_index,previous)
    vector=zeros(dimension,1);
    vector(mod(seed_index-1,dimension)+1)=1;
    vector=orthogonalize(vector,previous);
    if norm(vector)<=eps
        vector=orthogonalize((1:dimension)',previous);
    end
end

function vector=orthogonalize(vector,basis)
    if ~isempty(basis)
        vector=vector-basis*(basis'*vector);
        vector=vector-basis*(basis'*vector);
    end
    magnitude=norm(vector);
    if magnitude<=eps
        projector=eye(size(basis,1))-basis*basis';
        [~,index]=max(vecnorm(projector,2,1));
        vector=projector(:,index);
        magnitude=norm(vector);
    end
    vector=vector/max(magnitude,eps);
end

function [vector,value]=leading_eigenvector(matrix)
    matrix=(matrix+matrix')/2;
    [vectors,values]=eig(matrix,'vector');
    [value,index]=max(real(values));
    vector=real(vectors(:,index));
end

function vector=deterministic_sign(vector)
    [~,index]=max(abs(vector));
    if vector(index)<0, vector=-vector; end
end

function value=ssqcor_objective(cross,vectors)
    value=0;
    for m=1:numel(vectors)-1
        for n=m+1:numel(vectors)
            correlation=vectors{m}'*cross{m,n}*vectors{n};
            value=value+abs(correlation)^2;
        end
    end
end

function [weights,scores]=sort_mcca_components(views,weights)
    rank_value=size(weights{1},2);
    scores=zeros(rank_value,1);
    for k=1:rank_value
        component_sets=cellfun(@(x,w) x*w(:,k),views,weights,'UniformOutput',false);
        correlations=[];
        for m=1:numel(views)-1
            for n=m+1:numel(views)
                correlations(end+1,1)=absolute_correlation( ...
                    component_sets{m},component_sets{n}); %#ok<AGROW>
            end
        end
        scores(k)=sqrt(mean(correlations.^2));
    end
    [scores,order]=sort(scores,'descend');
    for m=1:numel(weights), weights{m}=weights{m}(:,order); end
end

function features=extract_component_features( ...
        scores,mixing,predictability,emg,Fs,peripheral_mask)
    num_components=size(scores,2);
    high_frequency_fraction=zeros(num_components,1);
    spectral_slope=zeros(num_components,1);
    emg_coherence=nan(num_components,1);
    peripheral_energy=zeros(num_components,1);
    spectral_window=min(size(scores,1),max(128,round(2*Fs)));
    spectral_window=max(32,spectral_window-mod(spectral_window,2));
    overlap=floor(spectral_window/2);
    nfft=max(512,2^nextpow2(spectral_window));
    for k=1:num_components
        [power,f]=pwelch(scores(:,k),hamming(spectral_window),overlap,nfft,Fs);
        total_mask=f>=1 & f<=min(100,Fs/2-1);
        high_mask=f>=30 & f<=min(100,Fs/2-1) & ~(f>=58 & f<=62);
        high_frequency_fraction(k)=sum(power(high_mask))/max(sum(power(total_mask)),eps);
        slope_mask=f>=20 & f<=min(100,Fs/2-1) & ~(f>=58 & f<=62) & power>0;
        if nnz(slope_mask)>=3
            fit=polyfit(log10(f(slope_mask)),log10(power(slope_mask)),1);
            spectral_slope(k)=fit(1);
        else
            spectral_slope(k)=NaN;
        end
        if ~isempty(emg)
            coherence=mscohere(scores(:,k),emg,hamming(spectral_window),overlap,nfft,Fs);
            emg_coherence(k)=mean(coherence(high_mask));
        end
        topography=mixing(k,:);
        peripheral_energy(k)=sum(abs(topography(peripheral_mask)).^2)/ ...
            max(sum(abs(topography).^2),eps);
    end
    features=table((1:num_components)',predictability,high_frequency_fraction, ...
        spectral_slope,emg_coherence,peripheral_energy, ...
        'VariableNames',{'Component','Predictability','HighFrequencyFraction', ...
        'SpectralSlope','EMGCoherence','PeripheralEnergy'});
end

function gains=automatic_component_gains(features,options)
    count=height(features);
    gains=ones(count,1);
    if options.MaxCandidateFraction==0 || count<2, return; end
    q_rank=percentile_rank(features.Predictability);
    hf_rank=percentile_rank(features.HighFrequencyFraction);
    slope_rank=percentile_rank(features.SpectralSlope);
    peripheral_rank=percentile_rank(features.PeripheralEnergy);
    if all(isnan(features.EMGCoherence))
        emg_rank=nan(count,1);
    else
        emg_rank=percentile_rank(features.EMGCoherence);
    end
    candidate_limit=max(1,floor(options.MaxCandidateFraction*count));
    [~,order]=sort(features.Predictability,'ascend');
    candidates=order(1:candidate_limit);
    for index=candidates(:)'
        votes=[hf_rank(index)>=0.75, slope_rank(index)>=0.75, ...
            peripheral_rank(index)>=0.75];
        if ~isnan(emg_rank(index))
            votes(end+1)=emg_rank(index)>=0.75 && ...
                features.EMGCoherence(index)>=0.05; %#ok<AGROW>
        end
        if sum(votes)<2, continue; end
        evidence=[1-q_rank(index),hf_rank(index),slope_rank(index), ...
            peripheral_rank(index)];
        if ~isnan(emg_rank(index)), evidence(end+1)=emg_rank(index); end %#ok<AGROW>
        confidence=mean(evidence);
        gains(index)=max(options.MinComponentGain, ...
            1-(1-options.MinComponentGain)*confidence^2);
    end
end

function ranks=percentile_rank(values)
    values=double(values(:));
    ranks=nan(size(values));
    valid=find(isfinite(values));
    if isempty(valid), return; end
    [sorted,order]=sort(values(valid));
    raw=zeros(numel(valid),1);
    first=1;
    while first<=numel(sorted)
        last=first;
        while last<numel(sorted) && sorted(last+1)==sorted(first), last=last+1; end
        raw(first:last)=mean([first last]);
        first=last+1;
    end
    if numel(valid)==1
        percent=zeros(size(raw));
    else
        percent=(raw-1)/(numel(valid)-1);
    end
    inverse=zeros(size(order)); inverse(order)=1:numel(order);
    ranks(valid)=percent(inverse);
end

function cleaned=reconstruct_with_gains(scores,mixing,residual,center,gains)
    cleaned=(scores.*gains(:)')*mixing+residual+center;
end

function metrics=preservation_metrics(original,cleaned,Fs)
    metrics=empty_preservation_metrics();
    if size(original,1)<max(64,round(2*Fs)), return; end
    upper=min(25,Fs/2-1);
    if upper<=0.5, return; end
    [b,a]=butter(4,[0.5 upper]/(Fs/2),'bandpass');
    before=filtfilt(b,a,original);
    after=filtfilt(b,a,cleaned);
    correlations=zeros(size(original,2),1);
    for c=1:size(original,2)
        correlations(c)=absolute_correlation(before(:,c),after(:,c));
    end
    before_power=median(mean(before.^2,1));
    after_power=median(mean(after.^2,1));
    covariance_before=cov(before);
    covariance_after=cov(after);
    triangle=triu(true(size(covariance_before)),1);
    if nnz(triangle)==0
        spatial=1;
    else
        spatial=absolute_correlation(covariance_before(triangle),covariance_after(triangle));
    end
    metrics.low_band_correlation=median(correlations);
    metrics.low_band_power_change_db=10*log10(max(after_power,eps)/max(before_power,eps));
    metrics.spatial_covariance_correlation=spatial;
end

function pass=passes_preservation(metrics,options)
    pass=isfinite(metrics.low_band_correlation) && ...
        metrics.low_band_correlation>=options.MinLowBandCorrelation && ...
        abs(metrics.low_band_power_change_db)<=options.MaxLowBandPowerChangeDb && ...
        metrics.spatial_covariance_correlation>=options.MinSpatialCovarianceCorrelation;
end

function metrics=empty_preservation_metrics()
    metrics=struct('low_band_correlation',NaN, ...
        'low_band_power_change_db',NaN, ...
        'spatial_covariance_correlation',NaN);
end

function diagnostics=empty_window_diagnostics()
    diagnostics=struct( ...
        'between_lag_scores',zeros(0,1), ...
        'component_features',table(), ...
        'component_gains',zeros(0,1), ...
        'component_unmixing',zeros(0,0), ...
        'component_mixing',zeros(0,0), ...
        'lag_unmixing',{{}}, ...
        'effective_eeg_rank',0, ...
        'mcca_objective',zeros(0,1), ...
        'mcca_sweeps',zeros(0,1), ...
        'reconstruction_error',0, ...
        'guard_strength',0, ...
        'preservation_before_guard',empty_preservation_metrics(), ...
        'preservation_after_guard',empty_preservation_metrics(), ...
        'skipped_reason',"");
end

function overrides=resolve_gain_overrides(requested,num_windows,num_channels)
    overrides=repmat({[]},num_windows,1);
    if isempty(requested), return; end
    if isnumeric(requested)
        validated=validate_gain_vector(requested,num_channels);
        overrides(:)={validated};
    elseif iscell(requested)
        if numel(requested)~=num_windows
            error('get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride', ...
                'A cell override must have one entry per window.');
        end
        for k=1:num_windows
            overrides{k}=validate_gain_vector(requested{k},num_channels);
        end
    else
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride', ...
            'ComponentGainsOverride must be numeric, a cell array, or empty.');
    end
end

function gains=validate_gain_vector(gains,num_channels)
    if ~isnumeric(gains) || ~isreal(gains) || ~isvector(gains) || ...
            any(~isfinite(gains)) || any(gains<0) || any(gains>1)
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride', ...
            'Component gains must be finite values in [0,1].');
    end
    gains=double(gains(:));
    if ~(isscalar(gains) || numel(gains)==num_channels)
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride', ...
            'A gain vector must be scalar or contain one gain per EEG channel.');
    end
end

function gains=validate_window_gain_override(requested,effective_rank)
    if isscalar(requested)
        gains=repmat(requested,effective_rank,1);
    elseif numel(requested)>=effective_rank
        gains=requested(1:effective_rank);
    else
        error('get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverrideForWindow', ...
            'The gain override is shorter than this window''s effective rank.');
    end
end

function metrics=calculate_metrics(original,cleaned,artifact,Fs,emg,channel_names)
    count=size(original,2);
    original_rms=sqrt(mean((original-mean(original,1)).^2,1))';
    cleaned_rms=sqrt(mean((cleaned-mean(cleaned,1)).^2,1))';
    artifact_rms=sqrt(mean(artifact.^2,1))';
    artifact_ratio=artifact_rms./max(original_rms,eps);
    preservation=zeros(count,1);
    low_band_correlation=zeros(count,1);
    low_band_power_change_db=zeros(count,1);
    high_frequency_power_change_db=zeros(count,1);
    emg_coherence_before=nan(count,1);
    emg_coherence_after=nan(count,1);
    [b,a]=butter(4,[0.5 min(25,Fs/2-1)]/(Fs/2),'bandpass');
    before_low=filtfilt(b,a,original);
    after_low=filtfilt(b,a,cleaned);
    spectral_window=min(size(original,1),max(128,round(2*Fs)));
    spectral_window=max(32,spectral_window-mod(spectral_window,2));
    overlap=floor(spectral_window/2);
    nfft=max(512,2^nextpow2(spectral_window));
    for c=1:count
        preservation(c)=absolute_correlation(original(:,c),cleaned(:,c));
        low_band_correlation(c)=absolute_correlation(before_low(:,c),after_low(:,c));
        low_band_power_change_db(c)=10*log10( ...
            max(mean(after_low(:,c).^2),eps)/max(mean(before_low(:,c).^2),eps));
        [p0,f]=pwelch(original(:,c),hamming(spectral_window),overlap,nfft,Fs);
        p1=pwelch(cleaned(:,c),hamming(spectral_window),overlap,nfft,Fs);
        high=f>=30 & f<=min(100,Fs/2-1) & ~(f>=58 & f<=62);
        high_frequency_power_change_db(c)=10*log10( ...
            max(sum(p1(high)),eps)/max(sum(p0(high)),eps));
        if ~isempty(emg)
            coh0=mscohere(original(:,c),emg,hamming(spectral_window),overlap,nfft,Fs);
            coh1=mscohere(cleaned(:,c),emg,hamming(spectral_window),overlap,nfft,Fs);
            emg_coherence_before(c)=mean(coh0(high));
            emg_coherence_after(c)=mean(coh1(high));
        end
    end
    metrics=table(channel_names,original_rms,cleaned_rms,artifact_rms,artifact_ratio, ...
        preservation,low_band_correlation,low_band_power_change_db, ...
        high_frequency_power_change_db,emg_coherence_before,emg_coherence_after, ...
        'VariableNames',{'Channel','OriginalRMS','CleanedRMS','ArtifactRMS', ...
        'ArtifactToSignalRMSRatio','OriginalToCleanedCorrelation','LowBandCorrelation', ...
        'LowBandPowerChangeDb','HighFrequencyPowerChangeDb', ...
        'EMGCoherenceBefore','EMGCoherenceAfter'});
end

function [matrix,channel_names]=eeg_structure_to_matrix(eeg_data)
    channel_names=fieldnames(eeg_data);
    if isempty(channel_names)
        error('get_eeg_windowed_mtlcca_artifact_removal:NoEEGChannels', ...
            'eeg_data must contain at least one channel.');
    end
    num_samples=numel(eeg_data.(channel_names{1}));
    matrix=zeros(num_samples,numel(channel_names));
    for c=1:numel(channel_names)
        signal=eeg_data.(channel_names{c});
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('get_eeg_windowed_mtlcca_artifact_removal:InvalidEEGChannel', ...
                'EEG channel "%s" must be a finite real numeric vector.',channel_names{c});
        end
        if numel(signal)~=num_samples
            error('get_eeg_windowed_mtlcca_artifact_removal:UnequalChannelLengths', ...
                'EEG channel "%s" has an inconsistent length.',channel_names{c});
        end
        matrix(:,c)=double(signal(:));
    end
end

function structure=matrix_to_eeg_structure(matrix,channel_names)
    structure=struct();
    for c=1:numel(channel_names)
        structure.(channel_names{c})=matrix(:,c);
    end
end

function starts=make_window_starts(num_samples,window_samples,hop_samples)
    if num_samples<=window_samples, starts=1; return; end
    last_start=num_samples-window_samples+1;
    starts=1:hop_samples:last_start;
    if starts(end)~=last_start, starts(end+1)=last_start; end
end

function weight=overlap_weight(num_samples)
    phase=pi*((1:num_samples)'-0.5)/num_samples;
    weight=sin(phase).^2;
end

function value=absolute_correlation(first,second)
    first=double(first(:))-mean(first);
    second=double(second(:))-mean(second);
    denominator=norm(first)*norm(second);
    if denominator<=eps
        if norm(first-second)<=eps(max([norm(first),norm(second),1]))
            value=1;
        else
            value=0;
        end
    else
        value=abs((first'*second)/denominator);
    end
end
