#include <gtest/gtest.h>
#include "blackbox_test.h"

int main(int argc, const char* argv[]) {
    set_running_path(argv[0]);
    ::testing::InitGoogleTest(&argc, (char**) argv);
    return RUN_ALL_TESTS();
}