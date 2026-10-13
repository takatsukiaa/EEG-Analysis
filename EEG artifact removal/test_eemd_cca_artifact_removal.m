function tests = test_eemd_cca_artifact_removal
%TEST_EEMD_CCA_ARTIFACT_REMOVAL Regression tests for EEMD-CCA.
%   results = runtests('test_eemd_cca_artifact_removal.m');
%   assertSuccess(results);
    tests = functiontests(localfunctions);
end

function testNoRejectionReconstructsSingleChannel(testCase)
    rng(201,'twister');
    Fs = 100;
    time = (0:599)'/Fs;
    eeg.Cz = sin(2*pi*8*time)+0.2*sin(2*pi*17*time)+0.03*randn(size(time));
    result = get_eemd_cca_artifact_removal(eeg,Fs, ...
        WindowDurationSec=6,OverlapFraction=0,EnsembleSize=2, ...
        MaxNumIMF=4,AutocorrelationThreshold=0,RandomSeed=11);
    verifyEqual(testCase,result.cleaned_eeg.Cz,eeg.Cz,'AbsTol',1e-9);
    verifyEqual(testCase,result.cleaned_eeg.Cz+result.artifact_estimate.Cz, ...
        eeg.Cz,'AbsTol',1e-10);
    verifyLessThan(testCase,max(result.eemd_reconstruction_error,[],'all'),1e-12);
    verifyLessThan(testCase,max(result.cca_reconstruction_error,[],'all'),1e-9);
end

function testEMGReferenceCorroboratesAutomaticSelection(testCase)
    rng(202,'twister');
    Fs = 200;
    time = (0:1599)'/Fs;
    emg = filter([1 -0.2],1,randn(size(time)));
    brain = sin(2*pi*8*time)+0.35*sin(2*pi*13*time+0.4);
    eeg.C3 = brain+1.2*emg;
    result = get_eemd_cca_artifact_removal(eeg,Fs, ...
        EMGReference=emg,WindowDurationSec=8,OverlapFraction=0, ...
        EnsembleSize=4,MaxNumIMF=6,AutocorrelationThreshold=0.95, ...
        EMGCorrelationThreshold=0.05,EMGCoherenceThreshold=0.05, ...
        MaxAutoRejectedFraction=0.5,RandomSeed=12);
    features = result.component_features{1,1};
    verifyTrue(testCase,any(features.EMGSupported));
    verifyNotEmpty(testCase,result.rejected_components{1,1});
    before = abs_correlation(eeg.C3,emg);
    after = abs_correlation(result.cleaned_eeg.C3,emg);
    verifyLessThan(testCase,after,before);
end

function testUnrelatedEMGPreventsAutomaticRejection(testCase)
    rng(203,'twister');
    Fs = 100;
    eeg.Fp1 = randn(700,1);
    unrelated = randn(700,1);
    result = get_eemd_cca_artifact_removal(eeg,Fs, ...
        EMGReference=unrelated,WindowDurationSec=7,OverlapFraction=0, ...
        EnsembleSize=2,MaxNumIMF=4,AutocorrelationThreshold=1, ...
        EMGCorrelationThreshold=1,EMGCoherenceThreshold=1, ...
        RequireEMGConfirmation=true,RandomSeed=13);
    verifyEmpty(testCase,result.rejected_components{1,1});
    verifyEqual(testCase,result.cleaned_eeg.Fp1,eeg.Fp1,'AbsTol',1e-9);
end

function testMultipleChannelsAndEMGStructure(testCase)
    rng(204,'twister');
    eeg.F3 = randn(600,1);
    eeg.F4 = randn(1,600);
    emg.Chin1 = randn(600,1);
    emg.Chin2 = randn(1,600);
    result = get_eemd_cca_artifact_removal(eeg,100, ...
        EMGReference=emg,WindowDurationSec=6,OverlapFraction=0, ...
        EnsembleSize=2,MaxNumIMF=3,AutocorrelationThreshold=0,RandomSeed=14);
    verifySize(testCase,result.cleaned_eeg.F3,[600 1]);
    verifySize(testCase,result.cleaned_eeg.F4,[600 1]);
    verifySize(testCase,result.component_features,[1 2]);
    verifyTrue(testCase,result.settings.has_emg_reference);
end

function testDeterministicWithoutChangingGlobalRandomStream(testCase)
    rng(205,'twister');
    eeg.C4 = randn(600,1);
    state_before = rng;
    first = get_eemd_cca_artifact_removal(eeg,100, ...
        WindowDurationSec=6,OverlapFraction=0,EnsembleSize=2, ...
        MaxNumIMF=4,AutocorrelationThreshold=0.9,RandomSeed=99);
    state_after = rng;
    second = get_eemd_cca_artifact_removal(eeg,100, ...
        WindowDurationSec=6,OverlapFraction=0,EnsembleSize=2, ...
        MaxNumIMF=4,AutocorrelationThreshold=0.9,RandomSeed=99);
    verifyEqual(testCase,state_after,state_before);
    verifyEqual(testCase,first.cleaned_eeg.C4,second.cleaned_eeg.C4,'AbsTol',0);
end

function testInvalidInputs(testCase)
    eeg.Cz = randn(200,1);
    verifyError(testCase,@() get_eemd_cca_artifact_removal( ...
        eeg,100,OverlapFraction=1), ...
        'get_eemd_cca_artifact_removal:InvalidOverlap');
    verifyError(testCase,@() get_eemd_cca_artifact_removal( ...
        eeg,100,SelectionDirection="sideways"), ...
        'get_eemd_cca_artifact_removal:InvalidSelectionDirection');
    verifyError(testCase,@() get_eemd_cca_artifact_removal( ...
        eeg,100,EMGReference=randn(199,1)), ...
        'get_eemd_cca_artifact_removal:InvalidEMGReference');
end

function value = abs_correlation(first,second)
    first = first(:)-mean(first);
    second = second(:)-mean(second);
    value = abs(first'*second)/(norm(first)*norm(second));
end
