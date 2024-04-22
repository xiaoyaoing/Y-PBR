#include "src/decode2d.c"

#include "constants/2dDouble.h"
#include "utils/rand64.h"
#include "zfpDecodeBlockBase.c"

int main()
{
  const struct CMUnitTest tests[] = {
    #include "testcases/block.c"
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
