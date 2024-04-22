#include "src/encode3l.c"

#include "constants/3dInt64.h"
#include "utils/rand64.h"
#include "zfpEncodeBlockBase.c"

int main()
{
  const struct CMUnitTest tests[] = {
    #include "testcases/block.c"
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
