function run_tests_v15_1_1()
%RUN_TESTS_V15_1_1 Entry point for the V15.1.1 prior-dimension hotfix.
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_1_1();
end
