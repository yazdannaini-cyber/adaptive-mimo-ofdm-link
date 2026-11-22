function run_tests_v15_1()
%RUN_TESTS_V15_1 Version-specific entry point, avoiding V15 path collisions.
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_1();
end
