//
// Sudoku Solver test suit
// by d/x, Spring/Summer/Autumn 2022, Daniel Koziarski
// TabSize = 2.
//

#include <cstdio>
#include "SudokuSolver.h"

#if 1
typedef struct {
  const char *name;
  uchar      clues[SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE];
} SudokuRiddle;

SudokuRiddle sudokuRiddles[] = {
  {
    "https://menneske.no  #0, Super easy, solution method: SiSo; 31 clues",
    {
      0,1,0, 0,9,3, 0,8,0,
      9,8,0, 0,0,0, 0,5,7,
      0,0,2, 0,0,8, 3,0,0,

      3,0,8, 9,0,0, 0,0,0,
      5,0,0, 0,1,0, 0,0,3,
      0,0,0, 0,0,2, 8,0,5,

      0,0,6, 4,0,0, 9,0,0,
      7,4,0, 0,0,0, 0,1,8,
      0,3,0, 5,8,0, 0,4,0, }
  },
  {
    "https://menneske.no  #1, Very easy, solution methods: SiSo, SC; 30 clues",
    {
      0,7,0, 9,0,8, 0,2,0,
      0,0,0, 3,0,6, 0,0,8,
      8,0,0, 0,5,0, 0,7,0,

      5,2,0, 0,0,0, 0,1,6,
      0,0,9, 0,0,0, 7,0,0,
      4,8,0, 0,0,0, 0,9,2,

      0,4,0, 0,2,0, 0,0,7,
      1,0,0, 8,0,7, 0,0,0,
      0,6,0, 4,0,3, 0,8,0, }
  },
  {
    "https://menneske.no  #2, Easy, solution methods: SiSo, SC, SB; 32 clues",
    {
      1,0,9, 0,0,7, 6,0,2,
      0,0,0, 6,0,1, 0,0,0,
      0,6,3, 0,9,0, 1,8,0,

      9,0,8, 0,0,0, 2,0,1,
      5,2,0, 0,0,0, 0,6,8,
      0,0,0, 0,0,0, 0,0,0,

      0,4,1, 0,2,0, 8,7,0,
      0,0,0, 9,0,6, 0,0,0,
      2,0,5, 7,0,0, 3,0,6, }
  },
  {
    "https://menneske.no  #3, Medium, solution methods: SiSo, SC, SB, DS; 26 clues",
    {
      3,0,6, 8,0,0, 2,0,9,
      0,0,0, 0,0,0, 0,0,0,
      5,0,0, 1,0,9, 0,0,4,

      0,0,7, 0,0,8, 1,0,6,
      0,0,0, 0,0,0, 0,0,0,
      2,0,4, 6,0,0, 3,0,0,

      6,0,0, 3,0,5, 0,0,7,
      0,0,0, 0,0,0, 0,0,0,
      4,0,9, 0,0,2, 8,0,3, }
  },
  {
    "https://menneske.no  #4, Hard, solution methods: SiSo, SC, SB, DS, DC; 23 clues",
    {
      0,0,0, 0,0,0, 0,0,0,
      4,0,8, 1,0,0, 0,0,7,
      0,3,0, 0,0,0, 0,0,0,

      0,0,0, 0,0,2, 8,0,0,
      6,0,1, 0,5,9, 0,0,0,
      0,0,7, 0,4,0, 9,0,0,

      0,0,0, 8,0,0, 0,9,0,
      0,0,0, 0,6,5, 0,2,0,
      0,0,5, 0,9,0, 0,8,6, }
  },
  {
    "https://menneske.no  #5, Harder, solution methods: SiSo, SC, SB, XW; 26 clues",
    {
      4,0,0, 0,0,0, 6,7,0,
      0,0,0, 8,0,0, 0,0,0,
      0,2,0, 0,0,4, 0,0,0,

      0,0,0, 7,9,0, 4,0,0,
      0,0,4, 0,3,6, 0,0,0,
      3,5,0, 4,0,0, 9,0,7,

      0,9,1, 0,0,0, 0,0,5,
      0,0,0, 0,0,2, 0,9,0,
      0,0,2, 5,0,0, 0,4,3, }
  },
  {
    "https://menneske.no  #6, Harder, solution methods: SiSo, SC, SB, DS, XW; 27 clues",
    {
      6,0,3, 0,0,2, 0,0,0,
      0,2,0, 0,5,6, 0,7,0,
      0,8,9, 0,0,0, 0,0,0,

      0,0,0, 9,0,0, 2,0,0,
      2,9,0, 0,8,0, 0,3,7,
      0,0,4, 0,0,7, 0,0,0,

      0,0,0, 0,0,0, 3,0,6,
      0,5,0, 0,4,1, 8,9,0,
      0,0,0, 0,0,9, 0,0,4, }
  },
  {
    "https://menneske.no  #7, Harder, solution methods: SiSo, SC, SB, DS, XW (2x); 24 clues",
    {
      0,0,0, 0,1,0, 0,0,5,
      0,5,0, 0,3,8, 4,2,0,
      9,2,1, 0,0,4, 0,0,0,

      2,0,0, 0,0,0, 7,0,0,
      0,0,0, 0,0,0, 0,0,9,
      0,0,0, 0,0,3, 0,0,0,

      0,7,0, 6,0,0, 9,0,0,
      6,0,8, 0,0,0, 0,1,0,
      5,0,0, 0,0,2, 0,0,8, }
  },
  {
    "https://menneske.no  #8, Harder, solution methods: SiSo, SC, SB, XW; 28 clues",
    {
      0,0,0, 0,6,0, 0,0,0,
      7,0,0, 0,4,0, 0,2,8,
      0,9,0, 7,0,5, 6,0,0,

      0,0,9, 4,0,0, 8,0,0,
      1,4,0, 0,0,0, 0,7,6,
      0,0,7, 0,0,2, 4,0,0,

      0,0,4, 8,0,3, 0,9,0,
      0,8,0, 0,1,0, 0,0,7,
      0,2,0, 0,9,0, 0,0,0, }
  },
  {
    "https://menneske.no  #9, Very Hard, solution methods: SiSo, SC, SB, FI; 24 clues",
    {
      2,0,0, 0,7,0, 8,0,5,
      0,0,0, 0,0,9, 3,0,7,
      0,0,0, 0,0,0, 0,0,0,

      0,0,8, 0,0,6, 0,0,4,
      0,0,6, 0,1,0, 5,0,0,
      0,1,4, 5,0,0, 0,0,9,

      0,0,1, 0,2,0, 0,0,0,
      0,0,0, 0,0,8, 0,0,3,
      7,0,0, 0,0,0, 4,9,0, }
  },
  {
    "https://menneske.no  #10, Super hard, solution methods: ?; 30 clues",
    {
      0,0,0, 0,3,0, 0,0,0,
      0,1,7, 2,0,0, 0,8,3,
      0,0,0, 6,0,9, 0,0,1,

      0,0,6, 0,1,0, 8,3,0,
      3,0,0, 5,0,7, 0,0,9,
      0,9,1, 0,8,0, 5,0,0,

      2,0,0, 7,0,5, 0,0,0,
      6,8,0, 0,0,3, 9,7,0,
      0,0,0, 0,9,0, 0,0,0, }
  },
  {
    "https://menneske.no  #11, Impossile, solution methods: ?; 23 clues, naked quintet test",
    {
      5,1,0, 4,0,0, 0,0,0,
      0,8,0, 0,0,5, 0,0,2,
      0,0,0, 0,8,0, 7,0,0,

      0,9,0, 3,0,0, 0,0,1,
      7,0,0, 0,0,0, 0,0,3,
      0,0,0, 0,6,0, 0,2,0,

      0,0,5, 8,3,0, 0,0,0,
      1,0,0, 0,0,0, 0,7,0,
      0,0,0, 0,5,4, 0,6,0, }
  },
  #if 1
  {
    "https://???  #12, Wikipedia: unsorted bruteForce worst case: 69e6 recursive calls; 17 clues, naked super test",
    {
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,3, 0,8,5,
      0,0,1, 0,2,0, 0,0,0,

      0,0,0, 5,0,7, 0,0,0,
      0,0,4, 0,0,0, 1,0,0,
      0,9,0, 0,0,0, 0,0,0,

      5,0,0, 0,0,0, 0,7,3,
      0,0,2, 0,1,0, 0,0,0,
      0,0,0, 0,4,0, 0,0,9, }
  },
  #endif
  {
    "https://www.youtube.com/watch?v=OmdLzZ2-2aQ&t=6s  #13, [level VH or higher], solution methods: swordfish, sevenfish",
    {
      6,0,0, 0,3,9, 0,0,8,
      0,9,0, 5,0,0, 2,0,0,
      0,8,0, 0,0,4, 0,0,9,

      2,0,0, 0,0,3, 6,0,0,
      0,6,8, 7,0,2, 1,4,0,
      0,0,1, 4,0,0, 0,0,2,

      8,0,0, 9,0,0, 0,0,0,
      0,0,6, 0,0,7, 0,2,0,
      5,0,0, 0,2,0, 0,0,1, }
  },
  {
    "https://www.youtube.com/watch?v=OmdLzZ2-2aQ&t=6s  #14, [level VH or higher], solution methods: swordfish",
    {
      0,2,0, 1,0,0, 6,9,0,
      0,6,0, 2,0,0, 0,0,5,
      3,7,0, 6,9,5, 2,8,0,

      0,4,0, 0,0,1, 8,0,2,
      1,3,0, 9,2,8, 0,5,0,
      2,8,0, 4,0,0, 0,3,0,

      6,5,2, 8,4,0, 0,1,7,
      7,1,0, 0,0,0, 0,2,0,
      0,9,3, 7,1,2, 5,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=OmdLzZ2-2aQ&t=6s  #15, [level VH or higher], solution methods: swordfish",
    {
      7,1,9, 0,8,2, 0,0,3,
      4,0,0, 0,0,0, 8,0,9,
      5,8,6, 3,0,0, 7,0,0,

      8,0,0, 0,0,0, 1,7,0,
      9,0,0, 0,2,0, 0,0,8,
      0,5,7, 0,0,0, 0,0,6,

      0,9,0, 0,0,4, 3,8,7,
      6,7,0, 0,0,0, 0,0,1,
      3,0,0, 9,7,0, 6,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=8lvq3bKYHu8 (8:40) #16, [level ?], solution methods: jellyfish test",
    {
      0,2,0, 6,9,7, 5,0,0,
      3,9,5, 4,1,2, 0,7,0,
      0,6,0, 3,5,8, 0,2,9,

      2,0,9, 0,0,6, 0,5,0,
      5,4,0, 2,3,9, 0,6,7,
      6,0,3, 5,0,1, 9,0,2,

      0,5,6, 1,2,3, 0,9,0,
      9,3,0, 0,6,4, 2,0,5,
      0,0,2, 9,0,5, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=8lvq3bKYHu8 (6:00)  #17, [level ?], solution methods: DC, jellyfish test",
    {
      6,0,0, 0,0,5, 0,2,0,
      8,0,0, 0,2,6, 4,5,1,
      0,2,0, 0,8,1, 0,7,6,

      7,6,0, 0,4,0, 2,0,5,
      9,0,0, 6,7,2, 0,4,3,
      4,0,2, 0,5,0, 0,6,7,

      2,4,0, 5,6,0, 7,1,0,
      1,8,6, 2,0,7, 5,0,4,
      0,0,0, 0,1,4, 6,0,2, }
  },
  {
    "https://www.youtube.com/watch?v=RZPIjZbcyNU (5:00)  #18, [level], solution methods: remote pair test (don't kill too much!)",
    {
      0,6,5, 1,2,0, 0,0,7,
      2,0,0, 6,0,4, 1,0,0,
      0,0,1, 0,0,7, 0,2,6,

      7,3,4, 2,6,5, 9,1,8,
      5,1,9, 0,0,0, 2,6,4,
      6,2,8, 4,1,9, 5,7,3,

      8,4,2, 3,0,6, 7,0,1,
      0,0,0, 0,0,1, 6,4,2,
      1,0,6, 0,4,2, 3,8,0, }
  },
  {
    "https://www.youtube.com/watch?v=RZPIjZbcyNU (6:55)  #19, [level], solution methods: remote pair test (2 chains)",
    {
      0,1,0, 6,2,0, 0,8,0,
      6,2,0, 0,8,0, 5,9,0,
      0,0,0, 9,0,0, 6,0,2,

      0,5,1, 7,0,8, 0,0,0,
      3,7,6, 2,1,4, 9,5,8,
      0,0,0, 5,0,3, 7,0,0,

      1,0,7, 0,5,2, 8,6,9,
      5,6,8, 1,7,9, 2,4,3,
      0,9,0, 8,0,6, 1,7,5, }
  },
  {
    "https://www.youtube.com/watch?v=RZPIjZbcyNU (9:55)  #20, [level], solution methods: remote pair test (long chain, 2 kill zones)",
    {
      2,4,1, 8,7,3, 6,5,9,
      8,3,6, 0,9,5, 0,7,1,
      9,7,5, 1,0,6, 3,8,0,

      0,6,2, 3,8,4, 1,9,0,
      0,8,3, 7,1,9, 0,0,6,
      0,1,9, 5,6,2, 8,3,0,

      3,5,7, 6,0,1, 9,0,8,
      6,2,4, 9,5,8, 7,1,3,
      1,9,8, 0,3,7, 0,6,0, }
  },
  {
    "https://www.youtube.com/watch?v=RZPIjZbcyNU (10:30)  #21, [level], solution methods: remote pair test (two chains, one with 2 kill zones)",
    {
      0,7,2, 5,4,6, 0,8,3,
      0,0,5, 2,0,0, 0,7,6,
      0,0,0, 0,9,0, 4,5,2,

      0,0,4, 0,0,0, 5,3,7,
      0,5,7, 3,0,4, 8,9,1,
      3,1,0, 0,7,5, 2,6,4,

      0,0,6, 0,5,0, 3,1,0,
      5,0,0, 0,0,1, 7,2,0,
      7,0,1, 0,3,2, 6,4,5, }
  },
  {
    "https://www.youtube.com/watch?v=RZPIjZbcyNU (11:05)  #22, [level], solution methods: remote pair test (long chain, 3 cells in kill zones)",
    {
      0,0,6, 7,0,0, 5,0,0,
      4,7,0, 0,3,0, 6,1,0,
      0,9,2, 6,0,0, 7,0,4,

      2,3,9, 1,0,0, 8,5,7,
      0,0,1, 2,7,3, 9,4,6,
      6,4,7, 9,8,5, 3,2,1,

      0,0,4, 3,0,7, 2,9,0,
      7,0,0, 4,9,0, 1,6,3,
      9,0,3, 0,0,0, 4,7,0, }
  },
  {
    "https://www.youtube.com/watch?v=krhYZS7uG3M (7:40)  #23, [level], solution methods: remote pair test (branchy chain)",
    {
      0,1,6, 8,4,5, 0,2,9,
      0,0,5, 7,9,3, 4,6,1,
      4,9,0, 1,6,2, 5,8,0,

      6,0,0, 9,3,4, 1,5,8,
      5,3,8, 2,1,6, 9,7,4,
      1,4,9, 5,7,8, 2,3,6,

      9,5,0, 6,8,1, 0,4,0,
      0,0,4, 3,2,9, 6,1,5,
      0,6,1, 4,5,7, 8,9,0, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (10:39) #24, [level], solution methods: empty rectangle test",
    {
      1,0,2, 0,4,6, 0,0,8,
      3,0,6, 8,0,0, 4,0,0,
      8,7,4, 1,0,0, 5,0,6,

      6,3,0, 0,0,0, 0,8,4,
      2,8,5, 4,0,9, 6,0,0,
      4,1,0, 0,6,8, 0,0,5,

      7,4,3, 0,0,5, 8,6,0,
      9,6,1, 0,8,4, 0,5,7,
      5,2,8, 6,0,0, 0,4,0, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (11:50) #25, [level], solution methods: empty rectangle test",
    {
      6,1,7, 8,0,0, 9,3,4,
      0,0,2, 3,6,0, 0,1,0,
      0,5,0, 1,0,4, 2,0,6,

      5,0,0, 0,1,8, 0,6,0,
      1,0,6, 4,0,0, 0,0,0,
      0,9,0, 6,3,0, 0,0,1,

      0,0,9, 2,4,1, 6,5,3,
      0,0,1, 5,8,6, 7,0,0,
      2,6,5, 0,0,3, 1,4,8, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (12:15) #26, [level], solution methods: empty rectangle test",
    {
      5,0,4, 9,0,0, 0,8,0,
      8,0,0, 7,0,0, 9,0,0,
      0,9,0, 0,8,1, 5,0,0,

      0,0,9, 0,0,0, 2,0,8,
      0,8,5, 0,0,0, 0,9,0,
      7,0,1, 8,0,9, 6,0,5,

      9,0,0, 1,4,0, 8,5,0,
      1,5,8, 0,0,3, 0,0,9,
      4,6,0, 5,9,8, 1,0,2, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (12:30) #27, [level], solution methods: empty rectangle test",
    {
      0,0,0, 0,0,0, 0,0,8,
      3,7,6, 0,0,0, 0,0,0,
      0,8,0, 0,2,1, 7,6,3,

      0,0,0, 3,6,0, 0,0,4,
      0,0,3, 0,0,0, 6,0,0,
      2,6,0, 0,4,9, 0,3,0,

      0,9,4, 5,1,8, 3,2,0,
      0,0,0, 0,0,0, 9,8,5,
      8,0,0, 0,9,0, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (12:40) #28, [level], solution methods: empty rectangle test",
    {
      0,7,0, 0,5,0, 6,8,4,
      5,1,8, 4,0,6, 2,0,0,
      0,4,6, 7,8,0, 0,1,5,

      6,0,7, 0,0,9, 5,4,2,
      0,0,5, 2,6,4, 1,0,0,
      4,2,1, 5,7,0, 0,0,6,

      1,6,0, 0,2,7, 4,5,0,
      0,0,0, 0,4,5, 7,6,1,
      7,5,4, 6,0,0, 0,2,0, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (12:50) #29, [level], solution methods: empty rectangle test",
    {
      0,8,0, 9,0,0, 1,7,6,
      7,1,2, 8,5,6, 3,9,4,
      6,0,0, 1,7,0, 2,8,5,

      0,0,8, 0,1,9, 0,0,0,
      2,0,6, 4,8,7, 0,1,3,
      1,0,0, 2,6,0, 8,0,0,

      0,0,3, 0,9,1, 0,0,8,
      9,0,0, 6,4,8, 0,3,0,
      8,6,0, 0,0,0, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=b2G8_0QE2ak&t=17s (13:15) #30, [level], solution methods: empty rectangle test",
    {
      1,0,2, 0,4,6, 0,0,8,
      3,0,6, 8,0,0, 4,0,0,
      8,7,4, 1,0,0, 5,0,6,

      6,3,0, 0,0,0, 0,8,4,
      2,8,5, 4,0,9, 6,0,0,
      4,1,0, 0,6,8, 0,0,5,

      7,4,3, 0,0,5, 8,6,0,
      9,6,1, 0,8,4, 0,5,7,
      5,2,8, 6,0,0, 0,4,0, }
  },
  {
    "https://www.youtube.com/watch?v=4bGzauHJmBM&t=9s (3:35) #31, [level], solution methods: unique rectangle typ1, naked super test",
    {
      1,0,0, 0,0,0, 0,0,4,
      0,0,5, 0,0,0, 2,0,1,
      0,2,0, 0,0,9, 0,6,0,

      0,0,1, 6,3,8, 7,0,0,
      8,5,7, 9,2,1, 4,3,6,
      2,3,6, 4,7,5, 9,1,8,

      0,8,0, 5,9,0, 1,4,0,
      7,0,9, 0,0,0, 5,0,0,
      5,0,0, 0,0,0, 6,0,9, }
  },
  {
    "https://www.youtube.com/watch?v=4bGzauHJmBM&t=9s (3:45) #32, [level], solution methods: unique rectangle typ1 test",
    {
      5,8,9, 7,4,0, 0,6,0,
      7,6,1, 9,2,0, 0,4,5,
      2,4,3, 6,5,0, 0,0,9,

      0,9,5, 3,8,7, 0,0,6,
      0,0,6, 0,1,9, 5,0,0,
      0,0,8, 0,6,5, 0,9,0,

      8,3,7, 1,9,2, 6,5,4,
      9,5,4, 8,3,6, 1,0,0,
      6,1,2, 5,7,4, 9,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=3r34Lh7DcuE&t=203s (3:05) #33, [level], solution methods: unique rectangle typ2 test",
    {
      3,0,9, 6,0,8, 5,4,7,
      8,5,0, 0,7,3, 0,0,2,
      6,7,0, 0,0,9, 8,3,1,

      1,0,3, 0,0,0, 0,5,9,
      7,0,0, 0,0,0, 1,0,3,
      2,9,0, 0,3,0, 4,0,6,

      4,0,7, 3,0,0, 0,0,8,
      9,3,0, 0,8,0, 0,1,5,
      5,0,0, 0,0,2, 3,0,4, }
  },
  {
    "https://www.youtube.com/watch?v=3r34Lh7DcuE&t=203s (4:05) #34, [level], solution methods: unique rectangle typ2 test",
    {
      0,0,1, 6,0,2, 0,7,0,
      0,5,0, 7,1,8, 0,0,0,
      2,0,0, 9,0,3, 0,1,6,

      1,6,8, 5,7,9, 0,0,0,
      5,0,0, 4,3,1, 6,8,7,
      0,0,0, 2,8,6, 1,5,9,

      9,0,0, 3,0,5, 0,4,1,
      0,1,0, 8,9,4, 0,6,0,
      0,4,0, 1,0,7, 3,9,0, }
  },
  {
    "https://www.youtube.com/watch?v=Yf15glek6yU&t=27s (5:20) #35, [level], solution methods: unique rectangle typ4 test",
    {
      1,0,4, 6,3,5, 0,0,7,
      0,6,0, 8,9,0, 0,0,0,
      0,0,0, 2,4,0, 0,6,0,

      6,5,8, 9,1,2, 7,3,4,
      3,0,0, 5,8,4, 6,1,9,
      4,1,9, 3,7,6, 2,8,5,

      0,0,0, 0,6,3, 0,7,2,
      0,0,0, 0,2,9, 0,5,0,
      2,0,0, 7,5,8, 1,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=Yf15glek6yU&t=27s (4:17) #36, [level], solution methods: unique rectangle typ4 test",
    {
      5,3,0, 0,0,7, 2,0,0,
      0,0,0, 0,2,0, 0,8,0,
      2,0,0, 4,0,3, 1,0,0,

      6,0,2, 0,0,4, 0,3,7,
      0,0,3, 0,7,0, 0,0,0,
      8,7,0, 3,0,0, 0,0,2,

      7,4,1, 6,3,9, 0,2,0,
      9,8,6, 0,4,0, 0,0,0,
      3,2,5, 7,0,0, 4,9,6, }
  },
  {
    "https://www.youtube.com/watch?v=jGFPrJeSY8Q&t=92s (4:10) #37, [level], solution methods: unique rectangle typ5 test",
    {
      5,1,0, 6,3,0, 0,7,0,
      0,6,3, 7,0,8, 0,0,0,
      0,0,7, 2,0,0, 6,0,3,

      0,0,9, 4,2,7, 0,3,6,
      3,0,0, 5,9,6, 0,0,7,
      7,0,6, 1,8,3, 4,9,0,

      0,7,0, 9,6,0, 3,0,0,
      0,3,0, 8,0,2, 7,6,0,
      6,0,0, 3,7,1, 0,2,4, }
  },
  {
    "https://www.youtube.com/watch?v=jGFPrJeSY8Q&t=92s (7:12) #38, [level], solution methods: unique rectangle typ5 test",
    {
      1,2,0, 3,7,0, 0,0,0,
      5,7,0, 2,8,1, 3,6,0,
      0,0,0, 0,4,5, 1,2,7,

      0,0,0, 0,0,4, 7,3,0,
      4,0,0, 0,0,0, 0,0,2,
      0,8,5, 0,0,3, 0,0,0,

      3,5,1, 4,0,8, 2,7,0,
      0,4,0, 0,3,0, 5,8,1,
      0,0,0, 1,5,2, 0,4,3, }
  },
  {
    "https://www.youtube.com/watch?v=jGFPrJeSY8Q&t=92s (5:00) #39, [level], solution methods: unique rectangle typ5 test",
    {
      1,6,0, 0,0,0, 4,7,0,
      8,0,3, 4,0,6, 2,0,0,
      0,0,5, 0,1,0, 6,0,3,

      0,2,0, 6,0,0, 0,0,0,
      9,3,1, 2,8,0, 0,6,4,
      0,8,6, 0,0,9, 0,2,0,

      6,0,8, 0,3,4, 9,0,2,
      0,0,0, 0,6,8, 0,3,7,
      3,5,0, 0,0,0, 8,4,6, }
  },
  {
    "https://www.youtube.com/watch?v=jGFPrJeSY8Q&t=92s (8:12) #40, [level], solution methods: unique rectangle typ5 test",
    {
      0,4,8, 5,1,2, 0,0,0,
      0,7,0, 9,3,8, 2,0,4,
      0,3,2, 4,6,7, 8,0,1,

      0,0,7, 1,8,9, 0,4,5,
      0,1,0, 2,5,3, 0,7,8,
      8,5,0, 6,7,4, 1,0,2,

      7,8,5, 3,9,1, 4,2,6,
      0,0,0, 8,2,5, 0,1,0,
      0,0,0, 7,4,6, 5,8,0, }
  },
  {
    "https://www.youtube.com/watch?v=5rqn0rR1x4A (7:35) #41, [level], solution methods: w-Wing test",
    {
      3,8,6, 0,4,0, 9,0,7,
      7,2,4, 8,0,9, 0,0,0,
      5,9,1, 7,3,0, 0,0,8,

      2,5,8, 4,7,1, 3,6,9,
      0,0,0, 0,8,0, 0,0,2,
      6,4,7, 9,2,3, 8,5,1,

      0,0,0, 0,5,0, 2,8,6,
      0,0,5, 2,0,8, 0,9,0,
      8,0,2, 0,9,0, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=T5cJepcjhIU&t=590s (12:25) #42, [level], solution methods: hidden rectangle (unique rectangle typ7 test",
    {
      4,0,6, 1,0,0, 0,5,2,
      0,3,1, 2,4,5, 6,0,9,
      5,2,0, 0,0,6, 4,1,3,

      0,5,0, 0,0,4, 3,0,1,
      3,6,4, 5,8,1, 2,9,7,
      1,0,0, 0,0,0, 5,4,0,

      9,4,7, 6,0,0, 0,3,5,
      6,1,3, 0,5,7, 9,2,0,
      2,8,5, 0,0,9, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (2:35) #43, [level], solution methods: unique rectangle typ3 test",
    {
      9,7,6, 4,2,8, 1,5,3,
      0,0,0, 0,6,0, 7,0,0,
      0,0,0, 0,0,0, 2,8,6,

      0,0,7, 5,8,6, 9,3,1,
      6,8,5, 9,3,1, 4,2,7,
      1,9,3, 2,7,4, 8,6,5,

      3,1,2, 0,0,0, 6,0,0,
      0,6,0, 0,0,0, 5,1,2,
      0,5,9, 6,1,2, 3,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (3:33) #44, [level], solution methods: unique rectangle typ3 test",
    {
      0,0,0, 7,4,6, 2,3,1,
      0,0,0, 9,8,1, 0,0,0,
      0,0,0, 5,2,3, 0,0,0,

      8,0,0, 4,7,0, 1,0,0,
      2,7,1, 8,3,5, 6,9,4,
      0,0,4, 6,1,0, 0,0,7,

      0,0,0, 1,6,4, 0,0,0,
      0,0,0, 3,9,8, 0,0,0,
      6,9,3, 2,5,7, 4,1,8, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (8:34) #45, [level], solution methods: unique rectangle typ3 test",
    {
      0,0,0, 2,7,6, 0,0,0,
      2,0,0, 0,9,0, 1,0,7,
      9,0,0, 8,0,0, 0,6,2,

      6,3,2, 1,4,8, 7,5,9,
      5,0,0, 9,6,2, 0,0,1,
      1,9,8, 7,5,3, 0,0,6,

      8,2,0, 0,0,9, 0,7,5,
      0,0,9, 0,8,0, 0,0,3,
      0,0,0, 6,2,0, 0,0,0, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (8:54) #46, [level], solution methods: unique rectangle typX test",
    {
      0,5,0, 2,3,0, 9,0,0,
      0,0,0, 0,0,0, 0,0,2,
      2,0,0, 0,0,8, 4,5,0,

      4,6,2, 3,8,1, 5,7,9,
      7,8,0, 0,9,0, 2,0,1,
      0,1,0, 0,2,7, 8,4,0,

      0,4,1, 8,0,2, 0,0,0,
      3,2,0, 0,0,0, 0,0,0,
      0,0,8, 0,5,3, 0,2,0, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (10:39) #47, [level], solution methods: unique rectangle typX test",
    {
      9,0,5, 7,0,4, 0,0,3,
      0,0,3, 0,6,9, 4,0,0,
      0,1,4, 3,0,8, 9,0,0,

      0,0,8, 0,0,0, 0,0,1,
      0,7,9, 8,3,1, 0,5,0,
      1,0,6, 0,0,0, 3,0,0,

      0,9,2, 0,0,0, 0,3,0,
      0,0,7, 9,8,3, 0,0,0,
      0,0,1, 4,7,2, 5,0,9, }
  },
  {
    "https://www.youtube.com/watch?v=C8VhAT7M6Gc&t=336s (11:20) #48, [level], solution methods: unique rectangle typX test",
    {
      0,0,0, 0,0,0, 0,9,8,
      8,1,3, 9,2,4, 7,6,5,
      0,9,0, 6,0,0, 0,3,1,

      0,2,0, 0,0,6, 9,8,0,
      9,0,0, 0,4,0, 6,0,0,
      0,6,8, 0,9,0, 0,1,0,

      2,0,0, 0,0,9, 0,4,0,
      5,4,9, 7,3,1, 8,2,6,
      0,8,0, 4,0,2, 0,0,9, }
  },
  {
    "https://www.youtube.com/watch?v=boZqrkk6BqI (11:20) #49, [level], finned swordfish test",
    {
      5,3,0, 0,0,2, 6,7,0,
      7,9,6, 0,5,0, 4,2,0,
      0,0,2, 0,0,7, 0,5,9,

      1,2,0, 7,3,0, 0,6,0,
      0,7,0, 2,1,6, 0,0,0,
      0,0,0, 0,4,8, 2,1,7,

      2,0,0, 0,0,0, 7,0,0,
      0,0,9, 0,7,0, 5,0,2,
      0,8,7, 0,2,0, 0,9,6, }
  },
  {
    "https://www.youtube.com/watch?v=boZqrkk6BqI (12:02) #50, [level], solution methods: finned swordfish test",
    {
      0,0,2, 0,4,5, 9,0,0,
      0,7,0, 2,0,9, 1,0,5,
      9,0,5, 1,0,8, 2,0,0,

      0,0,7, 5,0,2, 4,3,9,
      0,2,0, 8,0,3, 5,1,0,
      5,3,0, 0,0,4, 8,2,0,

      7,0,0, 9,0,6, 3,0,2,
      2,0,0, 0,0,1, 7,9,0,
      0,9,0, 0,2,7, 6,5,1, }
  },
  {
    "https://imgur.com/2wtVYqi (0:00) #51, [level], solution methods: X-chain test (closed loop, 10 links with nodes)",
    {
      0,5,0, 2,3,0, 4,0,6,
      3,0,0, 6,0,0, 0,9,0,
      6,0,0, 0,0,9, 3,8,0,
    //6,0,0, 0,0,0, 0,0,0,

      5,7,6, 9,0,0, 1,4,0,
      4,0,0, 0,7,6, 0,0,0,
      0,3,8, 0,0,4, 7,6,0,

      2,1,5, 4,0,0, 6,0,0,
      0,6,0, 0,0,2, 0,0,4,
      7,0,0, 0,6,5, 0,2,0, }
  },
  {
    "https://www.youtube.com/watch?v=bxf8JoS8cPY&t=4s (6:50) #52, [level], solution methods: X-chain with nice loop",
    {
      8,4,2, 0,5,0, 7,3,0,
      5,1,3, 0,7,4, 0,8,2,
      0,0,6, 0,3,0, 5,1,4,

      1,0,4, 3,8,7, 2,9,0,
      0,0,0, 0,1,0, 3,4,7,
      3,2,7, 4,9,0, 8,0,1,

      2,0,9, 0,4,0, 1,0,0,
      0,8,1, 5,2,0, 0,0,0,
      0,3,5, 0,6,0, 0,2,0, }
  },
  {
    "Grouped.xChain.5.links.#2.png (:) #53, grouped node x-chain test",
    {
      0,9,0, 7,1,0, 2,4,6,
      0,6,0, 8,9,0, 5,0,3,
      0,0,0, 6,0,0, 0,0,0,

      6,0,1, 4,3,0, 0,9,0,
      9,0,0, 1,7,0, 0,0,5,
      0,7,0, 0,6,9, 1,0,0,

      0,0,6, 0,0,1, 0,0,0,
      7,0,9, 0,0,6, 0,2,0,
      0,2,3, 0,4,7, 0,0,0, }
  },
  {
    "Grouped.xChain.5.links.2.weak.links.turned.in.group.nodes.png (:) #54, grouped node x-chain test",
    {
      1,0,0, 0,0,0, 8,4,0,
      0,0,0, 9,0,0, 0,0,0,
      0,5,2, 1,4,8, 3,9,0,

      0,6,0, 7,0,0, 0,5,9,
      0,0,0, 0,5,0, 0,0,0,
      2,0,5, 0,9,1, 0,3,0,

      0,0,6, 4,1,3, 5,8,0,
      5,0,0, 0,0,9, 4,0,3,
      4,8,3, 0,0,0, 9,0,1, }
  },
  {
    "Grouped.xChain.5.links.PERPENDICULAR.group.nodes.WEAK.link.png (:) #55, grouped node x-chain test",
    {
      0,0,1, 4,0,0, 0,0,3,
      5,0,0, 0,1,6, 0,0,0,
      0,8,2, 0,0,0, 6,0,0,

      0,0,0, 0,0,0, 0,5,7,
      0,0,7, 1,6,5, 3,0,0,
      3,9,5, 0,0,0, 0,6,0,

      0,0,4, 0,0,0, 9,3,0,
      0,0,0, 8,4,0, 0,0,6,
      7,0,0, 0,0,2, 4,0,0, }
  },
  {
    "Grouped.xChain.5.links.png (:) #56, grouped node x-chain test",
    {
      0,9,0, 7,1,0, 2,4,6,
      0,6,0, 8,9,0, 5,0,3,
      0,0,0, 6,0,0, 0,0,0,

      0,0,1, 4,3,0, 0,9,0,
      9,0,0, 1,7,0, 0,0,5,
      0,7,0, 0,6,9, 1,0,0,

      0,0,0, 0,0,1, 0,0,0,
      7,0,9, 0,0,6, 0,2,0,
      0,2,3, 0,4,7, 0,0,0, }
  },
  #if 0
  {
    " (:) #, grouped node x-chain test",
    {
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,

      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,

      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, }
  },
  {
    "https://menneske.no  #?, [level], solution methods: SiSo, [SC, SB, DS, DC, XW, FI]",
    {
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,

      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,

      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, }
  },
  #endif
};
#endif

//#define TEST_ALL
//#define TEST_AA
//#define TEST_DA

int main( int argc, char *argv[] ) {

  printf("Sudoku solver test suit\n");

  size_t inError = 0, unsolved = 0, solved = 0;

  for (int i = 51 /* 51, 9 */; i < elementsof( sudokuRiddles ); i++) {
    SudokuGrid sudokuGrid { sudokuRiddles[i].name, &sudokuRiddles[i].clues };

    switch (sudokuGrid.sudokuSolve( i )) {
    case -1: inError  += 1; break;
    case  0: unsolved += 1; break;
    case  1: solved   += 1; break;
    default: assert( 0 );
    }

    //if (i == 52)
    //  break;  // break the test suit after...
  }
  printf( "\nSudoku riddles count: %zu   Results:  inError: %zu,  unsolved: %zu,  solved: %zu\n",
          inError + unsolved + solved, inError, unsolved, solved );
  return 0;
}  // -------------------------------------------------------------------------------------------


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
