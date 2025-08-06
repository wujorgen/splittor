#include <gtest/gtest.h>

TEST(DummyTest, AlwaysPass) {
    EXPECT_TRUE(true);
}

TEST(SmartyTest, AlwaysFail) {
    EXPECT_TRUE(false);
}