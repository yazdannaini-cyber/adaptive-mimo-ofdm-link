function run_tests_v15_2_1()
%RUN_TESTS_V15_2_1 Entry point for V15.2.1 soft-decision FEC tests.
projectRoot=fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_2_1();
end
