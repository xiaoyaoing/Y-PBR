TEST_F(ARRAY_DIMS_SCALAR_TEST_ITERS, given_partialBlocks_when_incrementIterator_then_positionTraversesCorrectly)
{
  // force partial block traversal
  EXPECT_NE(0u, arr.size_x() % BLOCK_SIDE_LEN);
  EXPECT_NE(0u, arr.size_y() % BLOCK_SIDE_LEN);

  size_t totalBlocksX = (arr.size_x() + 3) / 4;
  size_t totalBlocksY = (arr.size_y() + 3) / 4;
  size_t totalBlocks = totalBlocksX * totalBlocksY;

  iter = arr.begin();
  for (size_t count = 0; count < totalBlocks; count++) {
    // determine if block is complete or partial
    size_t distanceFromEnd = arr.size_x() - iter.i();
    size_t blockLenX = distanceFromEnd < BLOCK_SIDE_LEN ? distanceFromEnd : BLOCK_SIDE_LEN;

    distanceFromEnd = arr.size_y() - iter.j();
    size_t blockLenY = distanceFromEnd < BLOCK_SIDE_LEN ? distanceFromEnd : BLOCK_SIDE_LEN;

    // ensure entries lie in same block
    size_t blockStartIndexI = iter.i();
    size_t blockStartIndexJ = iter.j();

    for (size_t j = 0; j < blockLenY; j++) {
      for (size_t i = 0; i < blockLenX; i++) {
        EXPECT_EQ(blockStartIndexI + i, iter.i());
        EXPECT_EQ(blockStartIndexJ + j, iter.j());
        iter++;
      }
    }
  }

//  EXPECT_EQ(arr.end(), iter); // triggers googletest issue #742
  EXPECT_TRUE(arr.end() == iter);
}

// const iterators

TEST_F(ARRAY_DIMS_SCALAR_TEST_ITERS, given_partialBlocks_when_incrementConstIterator_then_positionTraversesCorrectly)
{
  // force partial block traversal
  EXPECT_NE(0u, arr.size_x() % BLOCK_SIDE_LEN);
  EXPECT_NE(0u, arr.size_y() % BLOCK_SIDE_LEN);

  size_t totalBlocksX = (arr.size_x() + 3) / 4;
  size_t totalBlocksY = (arr.size_y() + 3) / 4;
  size_t totalBlocks = totalBlocksX * totalBlocksY;

  citer = arr.cbegin();
  for (size_t count = 0; count < totalBlocks; count++) {
    // determine if block is complete or partial
    size_t distanceFromEnd = arr.size_x() - citer.i();
    size_t blockLenX = distanceFromEnd < BLOCK_SIDE_LEN ? distanceFromEnd : BLOCK_SIDE_LEN;

    distanceFromEnd = arr.size_y() - citer.j();
    size_t blockLenY = distanceFromEnd < BLOCK_SIDE_LEN ? distanceFromEnd : BLOCK_SIDE_LEN;

    // ensure entries lie in same block
    size_t blockStartIndexI = citer.i();
    size_t blockStartIndexJ = citer.j();

    for (size_t j = 0; j < blockLenY; j++) {
      for (size_t i = 0; i < blockLenX; i++) {
        EXPECT_EQ(blockStartIndexI + i, citer.i());
        EXPECT_EQ(blockStartIndexJ + j, citer.j());
        citer++;
      }
    }
  }

//  EXPECT_EQ(arr.cend(), citer); // triggers googletest issue #742
  EXPECT_TRUE(arr.cend() == citer);
}
