function run_tests_v15_3()
%RUN_TESTS_V15_3 Entry point for V15.3 regression and joint-controller tests.
projectRoot=fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_3();
end
