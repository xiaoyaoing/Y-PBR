#include "src/encode2i.c"

#include "constants/2dInt32.h"
#include "utils/rand32.h"
#include "zfpEncodeBlockStridedBase.c"

int main()
{
  const struct CMUnitTest tests[] = {
    #include "testcases/blockStrided.c"
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
