//
// Sudoku Solver - engine
// by d/x, Spring/Summer/Autumn 2022, Daniel Koziarski
// TabSize = 2.
//

//#include <cstddef>
#include <cassert>
#include <cstdio>
#include <cstdarg>
#include <array>
#include <algorithm>
#include <map>
#include <utility>
#include <tuple>

#include "SudokuSolver.h"

#if 0
int SudokuGrid::Cell::CandidateTraits::singleSolution( void ) const {
  // Return cell clue (value), if only single clue remains possible, else return zero
  int count = 0;
  int clue;
  for (int i = 0; i < elementsof( clues_ ); i++)
    if (clues_[i] != 0)
      count += 1, clue = i;

  return (count == 1) ? clue + 1 : 0;
}  // -------------------------------------------------------------------------------------------
#endif

#if 0
int SudokuGrid::Cell::CandidateTraits::possibleClueCount( void ) const {
  int count = 0;
  for (int i = 0; i < elementsof( flag_ ); i++)
    if (flag_[i] != 0)
      count += 1;

  return count;
}  // -------------------------------------------------------------------------------------------
#endif

SudokuGrid::Possibilities &SudokuGrid::Possibilities::operator += ( const SudokuGrid::Cell &cell ) {

  const auto *candidateClue = cell.candidate_.clues_;
  for (auto &item : this->possibilities_)
    if (*candidateClue++ != 0) {
      item.count_ += 1;
      item.cell_   = const_cast<SudokuGrid::Cell *> (&cell);
    }
  return *this;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::Cell::candidateCrossOut( const CandidateTraits &candidateTraits ) {
  bool found = false;
  assert( unsolved() );
  assert( candidateUnprotected() );

  for (int clueIndex = 0; clueIndex < elementsof( candidateTraits.clues_ ); clueIndex++)
    if (candidatePossible( clueIndex )  &&  candidateTraits.candidatePossible( clueIndex )) {
      candidateCrossOut( clueIndex );
      found = true;
    }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::Cell::candidateCrossOutInProtectedCell( int clueIndex ) {
  if (candidateProtected() || candidateImpossible( clueIndex ))
    return false;

  candidateCrossOut( clueIndex );
  printf( "candidateCrossOutInProtectedCell(clue:%d) @[row:%d, col:%d]\n",
          clueIndex + 1, rowNumber( this ), colNumber( this ) );
  return true;
} // ----------------------------------------------------------------------------------------

bool SudokuGrid::Cell::candidateCrossOutInProtectedCell( const CandidateTraits &candidateTraits ) {
  if (candidateProtected())
    return false;

  return candidateCrossOut( candidateTraits );
} // ----------------------------------------------------------------------------------------

bool SudokuGrid::Cell::candidateFitsInto( const Cell *cell ) const {
  // return true, if all this->cell's candidates (possible clue values) are possible for cell too.
  assert( this->unsolved()  &&  cell->unsolved());

  for (int i = 0; i < elementsof( cell->candidate_.clues_ ); i++)
    if (this->candidatePossible( i )  &&  cell->candidateImpossible( i ))
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------

int  SudokuGrid::Cell::singleSolution( void ) const {
  // Return cell clue (value: [1 ... SUDOKU_GRID_SIZE]) if only single clue remains possible, else return zero

  int count = 0;
  int clue;
  for (int i = 0; i < elementsof( candidate_.clues_ ); i++)
    if (candidatePossible( i ))
      count += 1, clue = i;

  return (count == 1) ? clue + 1 : 0;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::Cell::cellInit( uchar clue, int row, int col,
                                 std::array<Cell, SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE> *grid ) {
  assert( row >= 0  &&  row < SUDOKU_GRID_SIZE);
  assert( col >= 0  &&  col < SUDOKU_GRID_SIZE);
  setClue( clue );

#if 1  // 26.10.2022
  if (candidateProtected()) {
    candidateProtection( false );
    assert( candidateUnprotected() );
  }
#else
  candidate_.crossOutProtected_ = false;
#endif
  const uchar startFlag = (clue != 0) ? 0 : 1;
  for (auto &item : candidate_.clues_)
    item = startFlag;

  // calculate start address of row/column/box, where grid[row,col] cell is located
#if 1
  // avoid nasty int-arithmetic overflow warnings
  rowOfCellsStart_ = &(*grid)[0] + row * SUDOKU_GRID_SIZE;
  colOfCellsStart_ = &(*grid)[0] + col;
  boxOfCellsStart_ = &(*grid)[0] + row / SUDOKU_BOX_SIZE * (SUDOKU_BOX_SIZE * SUDOKU_GRID_SIZE) +
                                   col / SUDOKU_BOX_SIZE * SUDOKU_BOX_SIZE;
#else
  rowOfCellsStart_ = &(*grid)[row * SUDOKU_GRID_SIZE];
  colOfCellsStart_ = &(*grid)[col];
  boxOfCellsStart_ = &(*grid)[row / SUDOKU_BOX_SIZE * (SUDOKU_BOX_SIZE * SUDOKU_GRID_SIZE) +
                              col / SUDOKU_BOX_SIZE * SUDOKU_BOX_SIZE];
#endif
}  // -------------------------------------------------------------------------------------------

SudokuGrid::~SudokuGrid() {
  delete grid_;
  //printf( "SudokuGrid destructor.\n" );
}  // -------------------------------------------------------------------------------------------

//SudokuGrid::SudokuGrid() :
//  riddleName_ { "?" },
//  clues_ { nullptr },
//  grid_ { nullptr } {
//  printf( "SudokuGrid() DEFAULT constructor, unexpected behaviour !!!\n" );
//}  // -------------------------------------------------------------------------------------------

SudokuGrid::SudokuGrid( const char *riddleName, const uchar (*clues)[SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE] ) :
  riddleName_ { riddleName } ,
  clues_ { clues } {
  grid_ = new std::array<Cell, SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE>;
  //printf( "SudokuGrid( riddleName, clues ) constructor.\n" );
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::candidateHouseProtection( bool state, HouseIterator &houseIterator ) {

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    houseIterator[i].candidateProtection( state );

}  // -------------------------------------------------------------------------------------------

void SudokuGrid::candidateGroupProtection( bool state, HouseIterator &groupIterator ) {

  for (int i = 0; i < SUDOKU_BOX_SIZE; i++)
    groupIterator[i].candidateProtection( state );

}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::candidateGroupProtected( HouseIterator &groupIterator ) {
  // return true, if ALL cells in group node are protected from cross out
  assert( groupIterator.isLineType() );

  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++)
    if (groupIterator[i].candidateUnprotected())
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::candidateGroupUnprotected( HouseIterator &groupIterator ) {
  // return true, if ALL cells in group node are NOT protected from cross out
  assert( groupIterator.isLineType() );

  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++)
    if (groupIterator[i].candidateProtected())
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::candidateGroupUnprotected( HouseIterator &groupIterator, int clueIndex ) {
  // return true, if ALL relevant cells in group node are NOT protected from cross out
  //              relevant cells are only those with candidatePossible( clueIndex ) == true
  assert( groupIterator.isLineType() );

  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++)
    if (groupIterator[i].candidateProtected()  &&  groupIterator[i].candidatePossible( clueIndex ))
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::candidateCrossOutInProtectedHouse( HouseIterator &houseIterator, int clueIndex ) {
  size_t crossOutCount = 0;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell &cell = houseIterator[i];
    if (cell.solved()  ||  cell.candidateProtected()  ||  cell.candidateImpossible( clueIndex ))
      continue;

    printf( "%s[row:%d, col:%d]", (crossOutCount == 0) ? "candidateCrossOutInProtectedHouse() @" : ":",
            rowNumber( &cell ), colNumber( &cell ) );
    cell.candidateCrossOut( clueIndex );
    crossOutCount += 1;
  }
  if (crossOutCount != 0)
    printf( "::clue:%d (%zu of %d times)\n", clueIndex + 1, crossOutCount, SUDOKU_GRID_SIZE );

  return (crossOutCount != 0) ? true : false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::candidateCrossOutInProtectedHouse( HouseIterator &houseIterator, const Cell::CandidateTraits &candidateTraits) {
  size_t crossOutCount = 0;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell &cell = houseIterator[i];
    if (cell.solved() || cell.candidateProtected())
      continue;
    if (cell.candidateCrossOut( candidateTraits )) {
      printf ( "%s[row:%d, col:%d]", (crossOutCount == 0) ? "candidateCrossOutInProtectedHouse() @" : ":",
               rowNumber( &cell), colNumber( &cell ) );
      crossOutCount += 1;
    }
  }
  if (crossOutCount != 0) {
    auto trait = [&candidateTraits]( int idx ) -> char { return (candidateTraits.clues_[idx] != 0) ? '1' + idx : '-'; };
    printf( "::traits:[%c%c%c%c%c%c%c%c%c] (%zu of %d times)\n",
            trait(0), trait(1), trait(2), trait(3), trait(4), trait(5), trait(6), trait(7), trait(8),
            crossOutCount, SUDOKU_GRID_SIZE );
  }

  return (crossOutCount != 0) ? true : false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::candidateCrossOutInProtectedHouseIfInAllNeighboursSight( HouseIterator &houseIterator,
                                                                          int clueIndex,
                                                                          const Cell * const *neighbours,
                                                                          size_t neighboursQuantity ) {
  size_t crossOutCount = 0;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell &cell = houseIterator[i];
    if (cell.solved()  ||  cell.candidateProtected()  ||  cell.candidateImpossible( clueIndex ))
      continue;
    if (cell.inAllNeighboursSight( neighbours, neighboursQuantity )) {

      printf( "%s[row:%d, col:%d]", (crossOutCount == 0) ? "candidateCrossOutInProtectedHouse() @" : ":",
              rowNumber( &cell ), colNumber( &cell ) );
      cell.candidateCrossOut( clueIndex );
      crossOutCount += 1;
    }
  }
  if (crossOutCount != 0)
    printf( "::clue:%d (%zu of %d times)\n", clueIndex + 1, crossOutCount, SUDOKU_GRID_SIZE );

  return (crossOutCount != 0) ? true : false;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::candidateCrossOutInRowColBox( Cell *cell ) {
  // Postfix operators don't work well with polymorphic classes.
  // Classes: RowIterator, ColIterator, BoxIterator are derived form HouseIterator, which is an abstract class.
  // Hence, it is not easy possible to implement postfix increment/decrent operators in derived classes:
  // https://stackoverflow.com/questions/37353387/incremental-operator-overload-in-an-abstract-class-c
  // 
  // For this reason prefix operators must be used instead or
  // each class derived from HouseIterator must provide own postfix operators :
  //
  int clueIndex = cell->getClueIndex();
#if 1
  RowIterator rowOfCells { cell->rowOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    rowOfCells++->candidateCrossOut( clueIndex );                    // postfix operator provided by derived class

  ColIterator colOfCells { cell->colOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    colOfCells++->candidateCrossOut( clueIndex );

  BoxIterator boxOfCells { cell->boxOfCellsStart_ };
  for (int i = 0; i < SUDOKU_BOX_SIZE * SUDOKU_BOX_SIZE; i++)
    boxOfCells++->candidateCrossOut( clueIndex );
#else
  RowIterator rowOfCells { cell->rowOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    rowOfCells->candidateCrossOut( clueIndex ), ++rowOfCells;        // prefix operator used instead of postfix

  ColIterator colOfCells { cell->colOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    colOfCells->candidateCrossOut( clueIndex ), ++colOfCells;        // prefix operator used instead of postfix

  BoxIterator boxOfCells { cell->boxOfCellsStart_ };
  for (int i = 0; i < SUDOKU_BOX_SIZE * SUDOKU_BOX_SIZE; i++)
    boxOfCells->candidateCrossOut( clueIndex ), ++boxOfCells;        // prefix operator used instead of postfix
#endif
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::establishClueAndCrossItOutInRowColBox( Cell *cell, int clue ) {
  assert( clue != 0 );
  cell->setClue( clue );
  for (auto &item : cell->candidate_.clues_)
    item = 0;
  candidateCrossOutInRowColBox( cell );
}  // -------------------------------------------------------------------------------------------

int  SudokuGrid::candidateCountInHouse( HouseIterator &houseIterator, int clueIndex ) {
  int candidateCount = 0;

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++)
    candidateCount += (int) houseIterator[i].candidatePossible( clueIndex );

  return candidateCount;
}  // -------------------------------------------------------------------------------------------

int  SudokuGrid::candidateCountInGroup( HouseIterator &groupIterator, int clueIndex ) {
  int candidateCount = 0;
#if 0
  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
    const Cell *cell = &groupIterator[i];
    printf( "candidateCountInGroup() cell @[row:%d, col:%d]\n", rowNumber( cell ), colNumber( cell ) );
    candidateCount += (int) cell->candidatePossible( clueIndex );
  }
#else
  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++)
    candidateCount += (int) groupIterator[i].candidatePossible( clueIndex );
#endif
  return candidateCount;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::initGrid( void ) {
  assert( grid_ != nullptr );
  assert( grid_->size() == capacityof( *clues_ ) );

  int row = 0, col = 0;

  auto cellIt = grid_->begin();
  for (const auto &clue : *clues_) {
    cellIt++->cellInit( clue, row, col, grid_ );           // transfer clues, init logic and row/col/box pointers

    if (col == SUDOKU_GRID_SIZE - 1)
      row += 1, col  = 0;
    else
      col += 1;
  }
  for (auto &cell : *grid_)
    if (cell.solved())
      candidateCrossOutInRowColBox( &cell );

  printSudokuGrid( "After initGrid():\n" );
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::singleSolution( void ) {
  // Return false, if new clue cannot be found
  // This method checks if a single clue value remains possible in each grid cell

  bool found = false;
  for (auto &cell : *grid_)
    if (cell.unsolved()) {
      auto clue = cell.singleSolution();
      if (clue != 0)
        found = true, establishClueAndCrossItOutInRowColBox( &cell, clue );
    }
  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::singleCellAnalysis( HouseIterator &&houseIterator ) {           /* aka: Hidden Single */
  // return false, if neither new clue can be established
  // Algorithm level: basic. Description:
  //   YouTube: "dxSudoku #2 Hidden Single"
  // 
  // Same algorithm can be iterated row by row and then column by column.
  //
  // HouseIterator is the base class for a row_ or column_ or box_iterator class.
  // A group of cells (cluster) consists of always SUDOKU_GRID_SIZE cells.
  // Cells within the particular (current) group are traversed according
  // to an indexed access arithmetic provided by the derived class.
  bool found = false;
  Possibilities mergedSum = { 0 };

  for (int cluster = 0; cluster < SUDOKU_GRID_SIZE; cluster++, houseIterator.nextHouse())
    for (int i = 0; i < SUDOKU_GRID_SIZE; i++)             // traverse a group of cells
      houseIterator[i].mergeCellPossibilities( &mergedSum );

  for (int i = 0; i < elementsof( mergedSum.possibilities_ ); i++)
    if (mergedSum.possibilities_[i].count_ == 1)           // clue possible in a single cell only?
      found = true, establishClueAndCrossItOutInRowColBox( mergedSum.possibilities_[i].cell_, i + 1 );

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::singleCell( void ) {
  // Return false, if new clue cannot be found
  // This method checks for the presence of a clue value possible only once within a group of cells.
  // A cell group consists of cells in the same row (or same column or same box) in sudoku grid.
  bool found = false;

  if (singleCellAnalysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (singleCellAnalysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (singleCellAnalysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


void SudokuGrid::singleBoxAnalysisInit( HouseIterator &houseIterator,
                                        Possibilities  (*mergedSumData)[SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE] ) {

  for (int box = 0; box < elementsof( *mergedSumData ); box++) {
    (*mergedSumData)[box].reset();
    for (int i = 0; i < SUDOKU_BOX_SIZE; i++)                        // traverse cells in box portions along row/column
      houseIterator[box * elementsof( *mergedSumData ) + i].mergeCellPossibilities( &(*mergedSumData)[box] );
  }
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::singleBoxAnalysis( HouseIterator &houseIterator,
                                    Possibilities  (*mergedSumData)[SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE],
                                    Possibilities *(*mergedSum)    [SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE] ) {
  // return false, if neither clue possibility can be crossed out in given box of cells
  // Algorithm description:
  // if a particular clue is possible twice or triply within the current row AND all possible
  // locations are within the same grid box, then this clue value can be surely crossed out 
  // from remaining rows inside of this box.
  // Same principle is appliable to grid columns.

  bool found = false;

  singleBoxAnalysisInit( houseIterator, mergedSumData );
  candidateHouseProtection( true, houseIterator );

  for (int box = 0; box < elementsof( *mergedSum ); box++) {
    for (int clueIndex = 0; clueIndex < elementsof( (*mergedSum)[0]->possibilities_ ); clueIndex++) {
      if (0 == (*mergedSum)[0]->possibilities_[clueIndex].count_)
        continue;                                          // clueIndex impossible in current box
      if (0 != std::count_if( *mergedSum + 1, *mergedSum + elementsof( *mergedSum ),
                              [&]( auto it ) { return it->possibilities_[clueIndex].count_ != 0; } ))
        continue;                                          // clueIndex sadly also possible in other boxes

      BoxIterator boxOfCells { (*mergedSum)[0]->possibilities_[clueIndex].cell_->boxOfCellsStart_ };
      if (candidateCrossOutInProtectedHouse( boxOfCells, clueIndex ))
        found = true;
    }
    std::rotate( *mergedSum, *mergedSum + 1, *mergedSum + elementsof( *mergedSum ) );
  }
  candidateHouseProtection( false, houseIterator );

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::singleBox( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  
  Possibilities mergedSumData[SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE] {};
  Possibilities *mergedSum[elementsof( mergedSumData )] {};

  for (int i = 0; i < elementsof( mergedSum ); i++)
    mergedSum[i] = &mergedSumData[i];

  RowIterator rowOfCells { (*grid_).data() };                        // traverse through rows in (*grid_)
  for (int row = 0; row < SUDOKU_GRID_SIZE; row++, rowOfCells.nextHouse())
    if (singleBoxAnalysis( rowOfCells, &mergedSumData, &mergedSum ))
      found = true;

  ColIterator colOfCells { (*grid_).data() };                        // traverse through columns
  for (int col = 0; col < SUDOKU_GRID_SIZE; col++, colOfCells.nextHouse())
    if (singleBoxAnalysis( colOfCells, &mergedSumData, &mergedSum ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

// disjoinSubsetAnalysis() implements detecting of naked constellations only in full version (333, 4444, 55555, 666666, ...) 
bool SudokuGrid::disjoinSubsetAnalysis( HouseIterator &houseIterator ) {
  // Return false, if new clue cannot be found in given box of cells.
  // Algorithm level: intermediate. Description:
  // if a particular clues PAIR      is the only possibility in 2 cells of a row, or    (aka: naked pair)
  // if a particular clues TRIPLE    is the only possibility in 3 cells of a row, or    (aka: naked triple)
  // if a particular clues QUADRUPLE is the only possibility in 4 cells of a row, or    (aka: naked quad / quartet)
  // if a particular clues QUINTUPLE is the only possibility in 5 cells of a row, or    (aka: naked quintet)
  // if a particular clues SIXFOLD   is the only possibility in 6 cells of a row, or    (aka: naked sextet)
  // if a particular clues SEVENFOLD is the only possibility in 7 cells of a row,       (aka: naked septed)
  // then these clues can be surely crossed out from remaining cells of the same row.
  // 
  // Algorithm level: intermediate. Description:
  //   YouTube: "dxSudoku #4 Naked Pair"
  //   YouTube: "dxSudoku #7, #8 Naked Triple"
  //   YouTube: "dxSudoku #57 Naked Quad"
  // 
  // Same principle is appliable also to grid columns and grid boxes.
  bool found = false;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {                                 // through cells in row/col/box group
    Cell &keyCell = houseIterator[i];
    if (keyCell.solved())
      continue;

    int sameCandidateCellCount = 0;
    for (int k = 0; k < SUDOKU_GRID_SIZE; k++) {
      Cell &cell = houseIterator[k];
      if (cell.unsolved()  &&  cell.candidateIdentical( &keyCell ))
        sameCandidateCellCount += 1;
    }
    if (sameCandidateCellCount != keyCell.candidateCount())
      continue;
    for (int k = 0; k < SUDOKU_GRID_SIZE; k++) {
      Cell &cell = houseIterator[k];
      if (cell.solved()  ||  cell.candidateIdentical( &keyCell ))
        continue;

      if (cell.candidateCrossOut( keyCell.candidate_ ))
        found = true;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::disjoinSubset( void ) {         /* aka: naked pair, naked triple, naked quartet, naked ... */
  // return false, if new clue cannot be found
  bool found = false;

  RowIterator rowOfCells { (*grid_).data() };                        // traverse through rows in (*grid_)
  for (int row = 0; row < SUDOKU_GRID_SIZE; row++, rowOfCells.nextHouse())
    if (disjoinSubsetAnalysis( rowOfCells ))
      found = true;

  ColIterator colOfCells { (*grid_).data() };                        // traverse through columns
  for (int col = 0; col < SUDOKU_GRID_SIZE; col++, colOfCells.nextHouse())
    if (disjoinSubsetAnalysis( colOfCells ))
      found = true;

  BoxIterator boxOfCells { (*grid_).data() };                        // traverse through boxes
  for (int box = 0; box < SUDOKU_GRID_SIZE; box++, boxOfCells.nextHouse())
    if (disjoinSubsetAnalysis( boxOfCells ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

//  disjoinChainAnalysis() implements partially "XYZ formation (XY-XYZ-YZ)" method  (aka: naked triple)
//  however, a naked triple like: "XYZ-XYZ-YZ" works pretty good too, but is NOT covered by current implementation
bool SudokuGrid::disjoinChainAnalysis( HouseIterator &houseIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm description:
  // if a possible clue's chain of type: [head:XY, body:XYZ, tail:YZ] exists in the same row
  // (where X, Y, Z stand for different clue values), then these clue values (X, Y. Z)
  // can be surely crossed out from remaining row cells.
  // 
  // Same principle is appliable to grid columns and grid boxes.
  bool found = false;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {                                 // through cells in row/col/box group
    Cell &bodyCell = houseIterator[i];
    if (bodyCell.solved()  ||
        bodyCell.candidateCount() != 3)
      continue;                                                                // chain body must have 3 possibilities
  
    for (int k = 0; k < SUDOKU_GRID_SIZE - 1; k++) {
      Cell &headCell = houseIterator[k];
      if (headCell.solved()  ||
          headCell.candidateCount() != 2  ||  headCell.candidateMisfitsInto( &bodyCell ))
        continue;
      for (k += 1; k < SUDOKU_GRID_SIZE; k++) {
        Cell &tailCell = houseIterator[k];
        if (tailCell.solved()  ||
            tailCell.candidateCount() != 2  ||  tailCell.candidateMisfitsInto( &bodyCell ))
          continue;
        if (headCell.candidateIdentical( &tailCell ))
          continue;

        // desired chain [XY,XYZ,YZ] established!
        chainProtection( true,  headCell, bodyCell, tailCell );
        if (candidateCrossOutInProtectedHouse( houseIterator, bodyCell.candidate_ ))
          found = true;
        chainProtection( false, headCell, bodyCell, tailCell );
        //if (found)
        //  return true;
      }
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::disjoinChain( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  RowIterator rowOfCells { (*grid_).data() };                        // traverse through rows in (*grid_)
  for (int row = 0; row < SUDOKU_GRID_SIZE; row++, rowOfCells.nextHouse())
    if (disjoinChainAnalysis( rowOfCells ))
      found = true;

  ColIterator colOfCells { (*grid_).data() };                        // traverse through columns
  for (int col = 0; col < SUDOKU_GRID_SIZE; col++, colOfCells.nextHouse())
    if (disjoinChainAnalysis( colOfCells ))
      found = true;

  BoxIterator boxOfCells { (*grid_).data() };                        // traverse through boxes
  for (int box = 0; box < SUDOKU_GRID_SIZE; box++, boxOfCells.nextHouse())
    if (disjoinChainAnalysis( boxOfCells ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

#if 0
bool SudokuGrid::nakedDuoAnalysis( HouseIterator &&houseIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: intermediate. Description:
  //   YouTube: "dxSudoku #4 Naked Pair"
  //
  //   Obviously possible type of naked duo: XY, XY.
  //   Same algorithm should be iterated row by row, then column by column, and then box by box.
  bool found = false;
  constexpr int N = 2;
  Cell *naked[N] {};

  for (int cluster = 0; cluster < SUDOKU_GRID_SIZE; cluster++, houseIterator.nextHouse())
    for (int i = 0; i < SUDOKU_GRID_SIZE - (N - 1); i++) {
      Cell *first = naked[0] = &houseIterator[i];
      if (first->solved() || first->candidateCount() > N)
        continue;
      for (int j = i + 1; j < SUDOKU_GRID_SIZE - (N - 2); j++) {
        Cell *second = naked[1] = &houseIterator[j];
        if (second->solved() || second->candidateCount() > N || second->candidateMisfitsInto( first ))
          continue;

        // naked set found
        chainProtection( true, naked, N );
        if (candidateCrossOutInProtectedHouse( houseIterator, /* killTraits */ second->candidate_ ))
          found = true;
        chainProtection( false, naked, N );
        break;
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------
#endif
bool SudokuGrid::nakedDuo( void ) {
  bool found = false;

  NakedMethodDuo nakedMethodDuo;

  if (nakedMethodDuo.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodDuo.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodDuo.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------
#if 0
bool SudokuGrid::nakedTrioAnalysis( HouseIterator &&houseIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: intermediate. Description:
  //   YouTube: "dxSudoku #7, #8 Naked Triple"
  // 
  //   Possible type of naked duo: XYZ-XYZ-XYZ, XYZ-XYZ-XY, XYZ-XZ-YZ, XY-XZ-YZ (all types in any combination)
  //   Same algorithm should be iterated row by row, then column by column, and then box by box.
  bool found = false;
  constexpr int N = 3;
  Cell *cell;
  Cell *naked[N] {};
  Cell::CandidateTraits joinedTraits[N];

  for (int cluster = 0; cluster < SUDOKU_GRID_SIZE; cluster++, houseIterator.nextHouse())
    for (int i = 0; i < SUDOKU_GRID_SIZE - (N - 1); i++) {
      cell = naked[0] = &houseIterator[i];
      if (cell->solved() || cell->candidateCount() > N)
        continue;
      joinedTraits[0] = cell->candidate_;

      for (int j = i + 1; j < SUDOKU_GRID_SIZE - (N - 2); j++) {
        cell = naked[1] = &houseIterator[j];
        if (cell->solved() || cell->candidateCount() > N)
          continue;
        if ((joinedTraits[1] = cell->candidate_ + joinedTraits[0]).candidateCount() > N)
          continue;

        for (int k = j + 1; k < SUDOKU_GRID_SIZE - (N - 3); k++) {
          cell = naked[2] = &houseIterator[k];
          if (cell->solved() || cell->candidateCount() > N)
            continue;
          if ((joinedTraits[2] = cell->candidate_ + joinedTraits[1]).candidateCount() > N)
            continue;

          // naked set found
          chainProtection( true, naked, N );
          if (candidateCrossOutInProtectedHouse( houseIterator, joinedTraits[2] ))
            found = true;
          chainProtection( false, naked, N );
          break;
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::nakedTrio( void ) {
  bool found = false;

  NakedMethodTrio nakedMethodTrio;

  if (nakedMethodTrio.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodTrio.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodTrio.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::nakedQuartet( void ) {
  bool found = false;
  NakedMethodQuartet nakedMethodQuartet;

  if (nakedMethodQuartet.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodQuartet.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodQuartet.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::nakedQuintet( void ) {
  bool found = false;
  NakedMethodQuintet nakedMethodQuintet;

  if (nakedMethodQuintet.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodQuintet.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodQuintet.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::nakedSextet( void ) {
  bool found = false;
  NakedMethodSextet nakedMethodSextet;

  if (nakedMethodSextet.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodSextet.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodSextet.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::nakedSeptet( void ) {
  bool found = false;
  NakedMethodSeptet nakedMethodSeptet;

  if (nakedMethodSeptet.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodSeptet.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodSeptet.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::nakedOctet( void ) {
  bool found = false;
  NakedMethodOctet nakedMethodOctet;

  if (nakedMethodOctet.analysis( RowIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodOctet.analysis( ColIterator { (*grid_).data() } ))
    found = true;
  if (nakedMethodOctet.analysis( BoxIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

SudokuGrid::Cell *SudokuGrid::lookBackForCellWithSamePotentialClue( HouseIterator &lineIterator, int clueIndex ) {

  while ((--lineIterator)->candidate_.clues_[clueIndex] == 0);

  return lineIterator.cell_;
}  // -------------------------------------------------------------------------------------------
#if 1
SudokuGrid::Cell *SudokuGrid::lookBackForCellWithSamePotentialClue( HouseIterator &houseIterator,
                                                                    Cell *seedCell, int clueIndex ) {
  houseIterator.cell_ = seedCell;

  while ((--houseIterator)->candidate_.clues_[clueIndex] == 0);

  return houseIterator.cell_;
}  // -------------------------------------------------------------------------------------------
#endif
#if 1
template<size_t N>
SudokuGrid::CellFishSequence<N>::CellFishSequence( HouseIterator &houseIterator, bool rowMode, int clueIndex )
  : rowMode_ { rowMode }  {
  Cell **fish = fish_;
  size_t fishCapacity = N;

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (houseIterator[i].candidatePossible( clueIndex ) && fishCapacity-- > 0) 
      *fish++ = &houseIterator[i];
}  // ---------------------------------------------------------------------------------------

template<size_t N>
SudokuGrid::CellFishSequence<N>::CellFishSequence( HouseIterator &houseIterator, bool rowMode,
                                                   Possibilities &possibilities, int clueIndex ) : rowMode_ { rowMode } {

  if (N == 2)
    assert( possibilities.possibilities_[clueIndex].count_ == 2 );
  else
    assert( possibilities.possibilities_[clueIndex].count_ >= 2 );

  fish_[0] = possibilities.possibilities_[clueIndex].cell_;
  fish_[1] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[0], clueIndex );

  if (N > 2)       // CellTrio or more
    if (possibilities.possibilities_[clueIndex].count_ > 2)
      fish_[2] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[1], clueIndex );
    else
      fish_[2] = nullptr;

  if (N > 3)       // CellQuartet or more
    if (possibilities.possibilities_[clueIndex].count_ > 3)
      fish_[3] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[2], clueIndex );
    else
      fish_[3] = nullptr;

  if (N > 4)       // CellQuintet or more
    if (possibilities.possibilities_[clueIndex].count_ > 4)
      fish_[4] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[3], clueIndex );
    else
      fish_[4] = nullptr;

  if (N > 5)       // CellSextet or more
    if (possibilities.possibilities_[clueIndex].count_ > 5)
      fish_[5] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[4], clueIndex );
    else
      fish_[5] = nullptr;

  static_assert( N <= 6 );   // further CellFishSequence class instantiations not implemented
}  // -------------------------------------------------------------------------------------------

template<size_t N>
void SudokuGrid::CellFishSequence<N>::CellFishSequenceInit( bool rowMode, Possibilities &possibilities, int clueIndex ) {
  rowMode_ = rowMode;

  static_assert( N >= 2);
  assert( possibilities.possibilities_[clueIndex].count_ >= 2 );

  // make sure possibilities.possibilities_[clueIndex].count_ has been not invalidated (after candidate crossing out):
  assert( possibilities.possibilities_[clueIndex].count_ ==
          candidateCountInHouse( (rowMode)
            ? (HouseIterator &) RowIterator { possibilities.possibilities_[clueIndex].cell_->rowOfCellsStart_ }
            : (HouseIterator &) ColIterator { possibilities.possibilities_[clueIndex].cell_->colOfCellsStart_ }, clueIndex ) );

  fish_[0] = possibilities.possibilities_[clueIndex].cell_;

  HouseIterator &lineIterator = (rowMode) ? (HouseIterator &) RowIterator { fish_[0] } 
                                          : (HouseIterator &) ColIterator { fish_[0] };
#if 1
  for (size_t i = 1; i < N; i++) {
    if (possibilities.possibilities_[clueIndex].count_ > i)
      fish_[i] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[i] = nullptr;
  }
#else
  fish_[1] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );

  if (N > 2)       // CellTrio or more
    if (possibilities.possibilities_[clueIndex].count_ > 2)
      fish_[2] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[2] = nullptr;

  if (N > 3)       // CellQuartet or more
    if (possibilities.possibilities_[clueIndex].count_ > 3)
      fish_[3] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[3] = nullptr;

  if (N > 4)       // CellQuintet or more
    if (possibilities.possibilities_[clueIndex].count_ > 4)
      fish_[4] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[4] = nullptr;

  if (N > 5)       // CellSextet or more
    if (possibilities.possibilities_[clueIndex].count_ > 5)
      fish_[5] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[5] = nullptr;

  if (N > 6)       // CellSeptet or more
    if (possibilities.possibilities_[clueIndex].count_ > 6)
      fish_[6] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[6] = nullptr;

  if (N > 7)       // CellOctet or more
    if (possibilities.possibilities_[clueIndex].count_ > 7)
      fish_[7] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    else
      fish_[7] = nullptr;
#endif

  static_assert( N <= 8 );   // further CellFishSequence class instantiations not implemented
}
#endif

//template class SudokuGrid::CellFishSequence<SudokuGrid::Cell, 2>;
#if 0
//how to implement methods in *.cpp for CellFishSequence template specialization without linking error?
template<>
SudokuGrid::CellSequence<SudokuGrid::Cell, 2>::CellSequence( HouseIterator &houseIterator, bool rowMode,
                                              Possibilities &possibilities, int clueIndex ) : rowMode_ { rowMode } {

  if (N == 2)
    assert( possibilities.possibilities_[clueIndex].count_ == 2 );
  else
    assert( possibilities.possibilities_[clueIndex].count_ >= 2 );

  fish_[0] = possibilities.possibilities_[clueIndex].cell_;
  fish_[1] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[0], clueIndex );

  if (N > 2)       // CellTrio or more
    if (possibilities.possibilities_[clueIndex].count_ > 2)
      fish_[2] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[1], clueIndex );
    else
      fish_[2] = nullptr;

  if (N > 3)       // CellQuartet or more
    if (possibilities.possibilities_[clueIndex].count_ > 3)
      fish_[3] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[2], clueIndex );
    else
      fish_[3] = nullptr;

  if (N > 4)       // CellQuintet or more
    if (possibilities.possibilities_[clueIndex].count_ > 4)
      fish_[4] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[3], clueIndex );
    else
      fish_[4] = nullptr;

  if (N > 5)       // CellSextet or more
    if (possibilities.possibilities_[clueIndex].count_ > 5)
      fish_[5] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[4], clueIndex );
    else
      fish_[5] = nullptr;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::x_Wing( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodDuo fishMethodDuo;

  if (fishMethodDuo.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodDuo.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::swordfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodTrio fishMethodTrio;

  if (fishMethodTrio.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodTrio.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

#if 1    // recursive templated version
bool SudokuGrid::jellyfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodQuartet fishMethodQuartet;

  if (fishMethodQuartet.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodQuartet.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------
#else    // without recursive template
bool SudokuGrid::jellyfishCrossOut( CellFishQuartet (&quartet)[4], int clueIndex ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  CellFishQuartet killZone = (*quartet)[0] + (*quartet)[1];
  for (size_t i = 2; i < capacityof( *quartet ); i++)
    killZone += (*quartet)[i];

  CellFishQuartet::chainProtection( true, quartet );
  for (size_t i = 0; i < elementsof( killZone.fish_ ); i++)
    if (killZone.fish_[i] != nullptr)
      if (killZone.rowMode_) {
        if (candidateCrossOutInProtectedHouse( ColIterator { killZone.fish_[i]->colOfCellsStart_}, clueIndex ))
          found = true;
      } else {
        if (candidateCrossOutInProtectedHouse( RowIterator { killZone.fish_[i]->rowOfCellsStart_}, clueIndex ))
          found = true;
      }
  CellFishQuartet::chainProtection( false, quartet );

  const bool rowMode = killZone.rowMode_;
  printf( "jellyfish(clues=%d) in %s: [%2d,%2d,%2d,%2d][%2d,%2d,%2d,%2d][%2d,%2d,%2d,%2d][%2d,%2d,%2d,%2d] kill zone [%s: %2d,%2d,%2d,%2d] (crossed out: %s)\n",
          clueIndex + 1, (rowMode) ? "rows" : "cols",
          ((*quartet)[0].fish_[0] != nullptr) ? (int) ((rowMode ? (*quartet)[0].fish_[0]->colOfCellsStart_ : (*quartet)[0].fish_[0]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[0].fish_[1] != nullptr) ? (int) ((rowMode ? (*quartet)[0].fish_[1]->colOfCellsStart_ : (*quartet)[0].fish_[1]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[0].fish_[2] != nullptr) ? (int) ((rowMode ? (*quartet)[0].fish_[2]->colOfCellsStart_ : (*quartet)[0].fish_[2]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[0].fish_[3] != nullptr) ? (int) ((rowMode ? (*quartet)[0].fish_[3]->colOfCellsStart_ : (*quartet)[0].fish_[3]->rowOfCellsStart_) - &(*grid_)[0]) : -1,

          ((*quartet)[1].fish_[0] != nullptr) ? (int) ((rowMode ? (*quartet)[1].fish_[0]->colOfCellsStart_ : (*quartet)[1].fish_[0]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[1].fish_[1] != nullptr) ? (int) ((rowMode ? (*quartet)[1].fish_[1]->colOfCellsStart_ : (*quartet)[1].fish_[1]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[1].fish_[2] != nullptr) ? (int) ((rowMode ? (*quartet)[1].fish_[2]->colOfCellsStart_ : (*quartet)[1].fish_[2]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[1].fish_[3] != nullptr) ? (int) ((rowMode ? (*quartet)[1].fish_[3]->colOfCellsStart_ : (*quartet)[1].fish_[3]->rowOfCellsStart_) - &(*grid_)[0]) : -1,

          ((*quartet)[2].fish_[0] != nullptr) ? (int) ((rowMode ? (*quartet)[2].fish_[0]->colOfCellsStart_ : (*quartet)[2].fish_[0]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[2].fish_[1] != nullptr) ? (int) ((rowMode ? (*quartet)[2].fish_[1]->colOfCellsStart_ : (*quartet)[2].fish_[1]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[2].fish_[2] != nullptr) ? (int) ((rowMode ? (*quartet)[2].fish_[2]->colOfCellsStart_ : (*quartet)[2].fish_[2]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[2].fish_[3] != nullptr) ? (int) ((rowMode ? (*quartet)[2].fish_[3]->colOfCellsStart_ : (*quartet)[2].fish_[3]->rowOfCellsStart_) - &(*grid_)[0]) : -1,

          ((*quartet)[3].fish_[0] != nullptr) ? (int) ((rowMode ? (*quartet)[3].fish_[0]->colOfCellsStart_ : (*quartet)[3].fish_[0]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[3].fish_[1] != nullptr) ? (int) ((rowMode ? (*quartet)[3].fish_[1]->colOfCellsStart_ : (*quartet)[3].fish_[1]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[3].fish_[2] != nullptr) ? (int) ((rowMode ? (*quartet)[3].fish_[2]->colOfCellsStart_ : (*quartet)[3].fish_[2]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          ((*quartet)[3].fish_[3] != nullptr) ? (int) ((rowMode ? (*quartet)[3].fish_[3]->colOfCellsStart_ : (*quartet)[3].fish_[3]->rowOfCellsStart_) - &(*grid_)[0]) : -1,

          (rowMode) ? "cols" : "rows",
          (killZone.fish_[0] != nullptr) ? (int) ((rowMode ? killZone.fish_[0]->colOfCellsStart_ : killZone.fish_[0]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          (killZone.fish_[1] != nullptr) ? (int) ((rowMode ? killZone.fish_[1]->colOfCellsStart_ : killZone.fish_[1]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          (killZone.fish_[2] != nullptr) ? (int) ((rowMode ? killZone.fish_[2]->colOfCellsStart_ : killZone.fish_[2]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          (killZone.fish_[3] != nullptr) ? (int) ((rowMode ? killZone.fish_[3]->colOfCellsStart_ : killZone.fish_[3]->rowOfCellsStart_) - &(*grid_)[0]) : -1,
          (found) ? "YES" : "-no-" );

  //if (found)
  //  printSudokuGrid( "after successful jellyfish cross out (clue: %d)\n", clueIndex + 1);

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::jellyfishAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   https://www.sudoku.org.pl/meduza.html
  //   YouTube: "dxSudoku #69 Jellyfish Puzzle Solving Technique"
  // Same principle is appliable to grid row and grid columns.
  //  ROW mode:
  //
  //  processing       fishHead[3] ---- fishHead[2] ---- fishHead[1] ---- fishHead[0]
  //  direction:           |                |                |
  //  top --> bottom   fishNeck[3] ---- fishNeck[2] ---- fishNeck[1] ---- fishNeck[0]
  //      row              |                |                |
  //                   fishBody[3] ---- fishBody[2] ---- fishBody[1] ---- fishBody[0]
  //                       |                |                |
  //                   fishTail[3] ---- fishTail[2] ---- fishTail[1] ---- fishTail[0]
  //
  //  Four jellyfish cell-quartets: ( fishHead[0,1,2,3], fishNeck[0,1,2,3], fishBody[0,1,2,3], fishTail[0,1,2,3] )
  //  However: each quartet have to contain at least 2 cells, to keep the jellyfish method working.
  //  A jellyfish row mode constellation consists of 8 to 16 cells, which occupy 4 rows and 2..4 columns.
  // 
  //  COLUMN mode:   (mode depends on houseIterator class)
  //
  //  processing       fishHead[3] ---- fishNeck[3] ---- fishBody[3] ---- fishTail[3]
  //  direction:           |                |                |
  //  left --> right   fishHead[2]      fishNeck[2] ---- fishBody[2]      fishTail[2]
  //      column           |                |                |
  //                   fishHead[1] ---- fishNeck[1] ---- fishBody[1] ---- fishTail[1]
  //                       |                |                |
  //                   fishHead[0] ---- fishNeck[0] ---- fishBody[0] ---- fishTail[0]
  //
  //  A jellyfish column mode constellation consists of 8 to 16 cells, which occupy 4 columns and 2..4 rows.
  //
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, houseIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      houseIterator[col].mergeCellPossibilities( &mergedSum[row] );

  CellFishQuartet quartet[4];

  CellFishQuartet &fishHead = quartet[0];
  CellFishQuartet &fishNeck = quartet[1];
  CellFishQuartet &fishBody = quartet[2];
  CellFishQuartet &fishTail = quartet[3];

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int headRow = 0; headRow < elementsof( mergedSum ) - 3; headRow++) {
      if (CellFishQuartet::invalidClueCount( mergedSum[headRow].possibilities_[clueIndex].count_ ))
        continue;
      fishHead.CellFishSequenceInit( houseIterator, rowMode, mergedSum[headRow], clueIndex );

      for (int neckRow = headRow + 1; neckRow < elementsof( mergedSum ) - 2; neckRow++) {
        if (CellFishQuartet::invalidClueCount( mergedSum[neckRow].possibilities_[clueIndex].count_ ))
          continue;
        fishNeck.CellFishSequenceInit( houseIterator, rowMode, mergedSum[neckRow], clueIndex );
        if (fishHead != fishNeck)                          // constellation invalid for jellyfish?
          continue;

        for (int bodyRow = neckRow + 1; bodyRow < elementsof( mergedSum ) - 1; bodyRow++) {
          if (CellFishQuartet::invalidClueCount( mergedSum[bodyRow].possibilities_[clueIndex].count_ ))
            continue;
          fishBody.CellFishSequenceInit( houseIterator, rowMode, mergedSum[bodyRow], clueIndex );
          if (fishHead + fishNeck != fishBody)             // constellation invalid for jellyfish?
            continue;

          for (int tailRow = bodyRow + 1; tailRow < elementsof( mergedSum ); tailRow++) {
            if (CellFishQuartet::invalidClueCount( mergedSum[tailRow].possibilities_[clueIndex].count_ ))
              continue;
            fishTail.CellFishSequenceInit( houseIterator, rowMode, mergedSum[tailRow], clueIndex );
            if (fishHead +fishNeck +fishBody != fishTail)  // constellation invalid for jellyfish?
              continue;

            printf( "jellyfish(clues=%d) clusters:%d,%d,%d,%d size:%d,%d,%d,%d (in %s)\n",
                    clueIndex + 1, headRow, neckRow, bodyRow, tailRow,
                    mergedSum[headRow].possibilities_[clueIndex].count_,
                    mergedSum[neckRow].possibilities_[clueIndex].count_,
                    mergedSum[bodyRow].possibilities_[clueIndex].count_,
                    mergedSum[tailRow].possibilities_[clueIndex].count_, (rowMode) ? "rows" : "cols" );
          #if 0
            if (FishMethodQuartet::fishCandidateCrossOut( quartet, killZone, clueIndex ))
              found = true;
          #else
            if (jellyfishCrossOut( quartet, clueIndex ))
              found = true;
          #endif
          }
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------
#endif

#if 1    // recursive templated version
bool SudokuGrid::starfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodQuintet fishMethodQuintet;

  if (fishMethodQuintet.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodQuintet.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------
#else    // without recursive template
bool SudokuGrid::starfishCrossOut( CellFishQuintet (&quintet)[5], int clueIndex ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  CellFishQuintet killZone = (*quintet)[0] + (*quintet)[1];
  for (size_t i = 2; i < capacityof( *quintet ); i++)
    killZone += (*quintet)[i];

  CellFishQuintet::chainProtection( true, quintet );
  for (size_t i = 0; i < elementsof( killZone.fish_ ); i++)
    if (killZone.fish_[i] != nullptr)
      if (killZone.rowMode_) {
        if (candidateCrossOutInProtectedHouse( ColIterator { killZone.fish_[i]->colOfCellsStart_}, clueIndex ))
          found = true;
      } else {
        if (candidateCrossOutInProtectedHouse( RowIterator { killZone.fish_[i]->rowOfCellsStart_}, clueIndex ))
          found = true;
      }
  CellFishQuintet::chainProtection( false, quintet );

  if (found)
    printf( "!!! starfish successful !!! (clue: %d)\n", clueIndex + 1 );

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::starfishAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   starfish goes one step further than jellyfish method (and two steps further than swordfish).
  //   An appliable YT video was not found.
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, houseIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      houseIterator[col].mergeCellPossibilities( &mergedSum[row] );

  CellFishQuintet killZone, quintet[5];          // N == 5

  CellFishQuintet &fishHead = quintet[0];
  CellFishQuintet &fishNeck = quintet[1];
  CellFishQuintet &fishBody = quintet[2];
  CellFishQuintet &fishLegs = quintet[3];
  CellFishQuintet &fishTail = quintet[4];

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int headRow = 0; headRow < elementsof( mergedSum ) - 4; headRow++) {
      if (CellFishQuintet::invalidClueCount( mergedSum[headRow].possibilities_[clueIndex].count_ ))
        continue;
      fishHead.CellFishSequenceInit( houseIterator, rowMode, mergedSum[headRow], clueIndex );

      for (int neckRow = headRow + 1; neckRow < elementsof( mergedSum ) - 3; neckRow++) {
        if (CellFishQuintet::invalidClueCount( mergedSum[neckRow].possibilities_[clueIndex].count_ ))
          continue;
        fishNeck.CellFishSequenceInit( houseIterator, rowMode, mergedSum[neckRow], clueIndex );
        if (fishHead != fishNeck)                          // constellation invalid for starfish?
          continue;
        killZone = fishHead + fishNeck;
        for (int bodyRow = neckRow + 1; bodyRow < elementsof( mergedSum ) - 2; bodyRow++) {
          if (CellFishQuintet::invalidClueCount( mergedSum[bodyRow].possibilities_[clueIndex].count_ ))
            continue;
          fishBody.CellFishSequenceInit( houseIterator, rowMode, mergedSum[bodyRow], clueIndex );
          if (killZone != fishBody)                        // constellation invalid for starfish?
            continue;
          killZone += fishBody;
          for (int legsRow = bodyRow + 1; legsRow < elementsof( mergedSum ) - 1; legsRow++) {
            if (CellFishQuintet::invalidClueCount( mergedSum[legsRow].possibilities_[clueIndex].count_ ))
              continue;
            fishLegs.CellFishSequenceInit( houseIterator, rowMode, mergedSum[legsRow], clueIndex );
            if (killZone != fishLegs)                      // constellation invalid for starfish?
              continue;
            killZone += fishLegs;
            for (int tailRow = legsRow + 1; tailRow < elementsof( mergedSum ); tailRow++) {
              if (CellFishQuintet::invalidClueCount( mergedSum[tailRow].possibilities_[clueIndex].count_ ))
                continue;
              fishTail.CellFishSequenceInit( houseIterator, rowMode, mergedSum[tailRow], clueIndex );
              if (killZone != fishTail)                    // constellation invalid for starfish?
                continue;
              killZone += fishTail;

              printf( "starfish(clues=%d) %s:%d,%d,%d,%d,%d (sizes:%d,%d,%d,%d,%d)\n",
                      clueIndex + 1, (rowMode) ? "rows" : "cols", headRow, neckRow, bodyRow, legsRow, tailRow,
                      mergedSum[headRow].possibilities_[clueIndex].count_,
                      mergedSum[neckRow].possibilities_[clueIndex].count_,
                      mergedSum[bodyRow].possibilities_[clueIndex].count_,
                      mergedSum[legsRow].possibilities_[clueIndex].count_,
                      mergedSum[tailRow].possibilities_[clueIndex].count_ );
            #if 0
              if (FishMethodQuintet::fishCandidateCrossOut( quintet, killZone, clueIndex ))
                found = true;
            #else
              if (starfishCrossOut( quintet, clueIndex ))
                found = true;
            #endif
            }
          }
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::starfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (starfishAnalysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (starfishAnalysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::sixfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodSextet fishMethodSextet;

  if (fishMethodSextet.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodSextet.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::sevenfish( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodSeptet fishMethodSeptet;

  if (fishMethodSeptet.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodSeptet.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::octopus( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;
  FishMethodOctet fishMethodOctet;

  if (fishMethodOctet.analysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (fishMethodOctet.analysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::w_WingCrossOut( const Cell *firstEP, const Cell *otherEP, int killCandidateIndex ) {
  bool found = false;
  const Cell *firstEndpointBox = firstEP->boxOfCellsStart_;
  const Cell *otherEndpointBox = otherEP->boxOfCellsStart_;

  if (firstEndpointBox->rowOfCellsStart_ == otherEndpointBox->rowOfCellsStart_  ||
      firstEndpointBox->colOfCellsStart_ == otherEndpointBox->colOfCellsStart_) {
    if (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( BoxIterator { firstEP->boxOfCellsStart_ },
                                                                 killCandidateIndex, &otherEP, 1))
      found = true;
    if (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( BoxIterator { otherEP->boxOfCellsStart_ },
                                                                 killCandidateIndex, &firstEP, 1))
      found = true;
  } else  // endpoints are NOT in same row or column of boxes, hence kill zone occurs in two single cells only
    if (intersectionOfCells( firstEP, otherEP )->candidateCrossOutInProtectedCell( killCandidateIndex )  ||
        intersectionOfCells( otherEP, firstEP )->candidateCrossOutInProtectedCell( killCandidateIndex ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::w_WingTraversing( RowIterator &&rowIterator, Cell *(*quartet)[4], int (*candidateXY)[2] ) {
  bool found = false;
  int  traversalIndex = (*candidateXY)[0];
  auto &[firstEP, otherEP, third] = *reinterpret_cast<const Cell *const(*)[3]> (*quartet);

  for (size_t row = 0; row < SUDOKU_GRID_SIZE; row++, rowIterator.nextHouse())
    for (size_t col = 0; col < SUDOKU_GRID_SIZE; col++) {
      const Cell *f_th /* = (*quartet)[3] */ = &rowIterator[col];
      if (f_th->solved() || f_th->candidateImpossible( traversalIndex ))
        continue;
      if (f_th == firstEP || f_th == otherEP || f_th == third)
        continue;
      if (f_th->occupyingDifferentHouse( third ) || f_th->occupyingDifferentHouse( otherEP ))
        continue;

      if (f_th->sharingSameRow( third ) && candidateCountInHouse( RowIterator { f_th->rowOfCellsStart_ }, traversalIndex ) > 2 ||
          f_th->sharingSameCol( third ) && candidateCountInHouse( ColIterator { f_th->colOfCellsStart_ }, traversalIndex ) > 2 ||
          f_th->sharingSameBox( third ) && candidateCountInHouse( BoxIterator { f_th->boxOfCellsStart_ }, traversalIndex ) > 2)
        continue;            // fourth and third are NOT an "either-or" link

      if (w_WingCrossOut( firstEP, otherEP, (*candidateXY)[1] )) {
        printf( "w-Wing  T-candidate:%d endpoints @[row:%d, col:%d][row:%d, col:%d]  T-path @[row:%d, col:%d][row:%d, col:%d]\n",
                traversalIndex + 1,
                rowNumber( firstEP ), colNumber( firstEP ), rowNumber( otherEP ), colNumber( otherEP ),
                rowNumber( third ),   colNumber( third ),   rowNumber( f_th ),    colNumber( f_th ) );
        found = true;
      }
    }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::w_WingTraversing( RowIterator &&rowIterator, Cell *(*quartet)[4] ) {
  bool  found         = false;
  int   candidateXY[2] {};
  const Cell *firstEP = (*quartet)[0];
//const Cell *otherEP = (*quartet)[1];

  for (int i = 0; i < elementsof( candidateXY ); i++)
    candidateXY[i] = firstEP->candidateNextIndex( (i == 0) ? -1 : candidateXY[i - 1] );

  for (size_t row = 0; row < SUDOKU_GRID_SIZE; row++, rowIterator.nextHouse())
    for (size_t col = 0; col < SUDOKU_GRID_SIZE; col++) {
      const Cell *third = (*quartet)[2] = &rowIterator[col];
      if (third->solved() || third == firstEP || third->occupyingDifferentHouse( firstEP ))
        continue;
      assert( third != (*quartet)[1] );
      for (int i = 0; i < elementsof( candidateXY ); i++) {          // try both candidates
        if (i > 0)
          std::rotate( candidateXY, candidateXY + 1, candidateXY + elementsof( candidateXY ) );

        if (third->candidateImpossible( candidateXY[0] ))
          continue;                                                  // does NOT match w-Wing sequence

        if (w_WingTraversing( RowIterator { third->rowOfCellsStart_->colOfCellsStart_ }, quartet, &candidateXY ))
          found = true;
      }
    }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::w_WingAnalysis( RowIterator &&rowIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #10 Revised W-Wing Puzzle Solving Technique"
  //   YouTube: "dxSudoku #10 Old Version"
  //
  // Analysis steps:
  // 1. Find first endpoint of w-Wing, cell type: XY (exact two candidates)
  // 2. Find other endpoint of w-Wing (same XY two candidates), which must be in a _different_ house (row/col/box).
  // 3. Find a cell (3rd) with X-candidate, which shares same house with the first endpoint.
  // 4. Find a cell (4th) with X-candidate, which shares same house with the "3rd" cell.
  //    The house shared by 3rd and 4th cell must contain only 2 cells with the X-candidate (-> either-or type link).
  // 5. The desired w-Wing occurs, if the 4th cell shares a house with the second endpoint (found in step 2).
  //    Traversing between the endpoints must have exactly 3 links (through 2 cells different from each other and
  //    from both endpoints).
  // 6. The Y-candidate can be crossed out from the kill zones, which are defined by cells shared by houses 
  //    of both w-Wing endpoints. 
  // 
  //   It makes no sense to iterate row by row _AND_THEN_ column by column when looking for w-Wing constellations.
  bool found = false;
  Cell *quartet[4] {};

  for (size_t row = 0; row < SUDOKU_GRID_SIZE - 1; row++, rowIterator.nextHouse()) 
    for (size_t col = 0; col < SUDOKU_GRID_SIZE; col++) {
      const Cell *firstEP = quartet[0] = &rowIterator[col];
      if (firstEP->solved() || firstEP->candidateCount() != 2)
        continue;
      RowIterator otherEndpointIterator { &rowIterator[0] };
      for (size_t otherRow = row + 1; otherRow < SUDOKU_GRID_SIZE; otherRow++) {
        otherEndpointIterator.nextHouse();
        for (size_t otherCol = 0; otherCol < SUDOKU_GRID_SIZE; otherCol++) {
          const Cell *otherEP = quartet[1] = &otherEndpointIterator[otherCol];
          if (otherEP->solved() || otherEP->candidateCount() != 2 || otherEP->candidateDifferent( firstEP ))
            continue;                                      // endpoints must provide same candidates of type XY
          if (otherEP->sharingSameBox( firstEP ) || otherEP->sharingSameCol( firstEP ))
            continue;                                      // endpoints must occupy different houses

          if (w_WingTraversing( RowIterator { firstEP->rowOfCellsStart_->colOfCellsStart_ }, &quartet ))
            found = true;
        }
      }
    }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::w_Wing( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (w_WingAnalysis( RowIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::xy_WingCrossOut( const Cell *pivot, const Cell *firstPincer, const Cell *otherPincer, int pincerZ ) {
  bool found = false;

  pivot->candidateProtection( true );
  if (w_WingCrossOut( firstPincer, otherPincer, pincerZ ))
    found = true;
  pivot->candidateProtection( false );

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xy_WingAnalysis( HouseIterator &otherPincerIterator,
                                  const Cell *pivot, const int (*pivotXY)[2], const Cell *firstPincer ) {
  bool found  = false;
  auto &[pivotX, pivotY] = *pivotXY;

  if (firstPincer->candidateImpossible( pivotX ) || firstPincer->candidatePossible( pivotY ))
    return false;                                          // firstPincer _NOT_ in  XZ formation

  auto otherPincerFormYZ = (firstPincer->candidate_ - pivotX) + pivotY;
  for (int k = 0; k < SUDOKU_GRID_SIZE; k++) {
    const Cell *otherPincer = &otherPincerIterator[k];
    if (otherPincer->solved() || otherPincer->candidateCount() != 2 || otherPincer == pivot)
      continue;
    if (otherPincer->candidateDifferent( otherPincerFormYZ ) || otherPincer->sharingSameHouse( firstPincer ))
      continue;

    int pincerZ = otherPincerFormYZ - pivotY;              // subtraction produces here only one candidate
    if (xy_WingCrossOut( pivot, firstPincer, otherPincer, pincerZ )) {
      printf( "XY-Wing: pivot(clues:%d,%d) @[row:%d, col:%d]  pincer(clue:%d) @[row:%d, col:%d]:[row:%d, col:%d]\n",
              pivotX + 1, pivotY + 1, rowNumber( pivot ), colNumber( pivot ),
              pincerZ + 1, rowNumber( firstPincer ), colNumber( firstPincer ),
              rowNumber( otherPincer ), colNumber( otherPincer ) );
      found = true;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xy_WingAnalysis( RowIterator &&rowIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #35 XY-Wing Puzzle Solving Technique"
  //
  // Analysis steps:
  // 1. Find pivot cell of type: XY (exact 2 candidates)
  //    The pivot cell has as usual three houses: RH,CH.BH (row, column, box)
  // 2. Find first pincer cell in form: XZ (exact 2 candidates) within first (of three) pivot houses
  // 3. Find other pincer cell in form: YZ (exact 2 candidates) within other (of remaining two) pivot houses,
  //    both pincer cells must be in different houses
  // 4. The Z-candidate can be crossed out from the kill zones, which are defined by every
  //    cell sharing a house with both xy-Wing pincers.
  // 
  //   It makes no sense to iterate row by row _AND_THEN_ column by column when looking for xy-Wing constellations.
  bool found = false;
  int  pivotXY[2] {};

  for (size_t row = 0; row < SUDOKU_GRID_SIZE; row++, rowIterator.nextHouse()) 
    for (size_t col = 0; col < SUDOKU_GRID_SIZE; col++) {
      const Cell *pivot = &rowIterator[col];
      if (pivot->solved() || pivot->candidateCount() != 2)
        continue;

      for (int i = 0; i < elementsof( pivotXY ); i++)
        pivotXY[i] = (i == 0) ? pivot->candidateFirstIndex() : pivot->candidateNextIndex( pivotXY[i - 1] );
      // works Ok, (casting required):
      HouseIterator *pincerIterators[3] = { &(RowIterator &) RowIterator { pivot->rowOfCellsStart_ },
                                            &(ColIterator &) ColIterator { pivot->colOfCellsStart_ },
                                            &(BoxIterator &) BoxIterator { pivot->boxOfCellsStart_ } };
      for (size_t i = 0; i < elementsof(pincerIterators); i++) {
        if (i > 0)
          std::rotate( pincerIterators, pincerIterators + 1, pincerIterators + elementsof( pincerIterators ) );

        auto &pincerIterator = *pincerIterators[0];
        for (int j = 0; j < SUDOKU_GRID_SIZE; j++) {
          const Cell *firstPincer = &pincerIterator[j];
          if (firstPincer->solved() || firstPincer->candidateCount() != 2 || firstPincer == pivot)
            continue;

          for (int r = 0; r < elementsof( pivotXY ); r++) {      // try both candidates
            if (r > 0)
              std::rotate( pivotXY, pivotXY + 1, pivotXY + elementsof( pivotXY ) );

            if (xy_WingAnalysis( *pincerIterators[1], pivot, &pivotXY, firstPincer ))
              found = true;
          }
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xy_Wing( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (xy_WingAnalysis( RowIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xyz_WingAnalysis( HouseIterator &&otherPincerIterator,
                                   Cell *(*trio)[3], Cell::CandidateTraits &otherPincerForm ) {
  bool found = false;
  auto &[pivot, firstPincer] = *reinterpret_cast<const Cell *const(*)[2]> (*trio);

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    const Cell *otherPincer = (*trio)[2] = &otherPincerIterator[i];
    if (otherPincer == pivot || otherPincer->sharingSameBox( pivot ) || otherPincer->solved())
      continue;
    if (otherPincer->candidateCount() != 2 || otherPincer->candidateDifferent( otherPincerForm ))
      continue;

    chainProtection( true, *trio, elementsof( *trio ) );
    int killClueIndex = pivot->candidate_ * firstPincer->candidate_ * otherPincer->candidate_;
    if (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( BoxIterator { pivot->boxOfCellsStart_ },
                                                                 killClueIndex, &(*trio)[1], elementsof( *trio ) - 1 )) {
      printf( "XYZ-Wing pivot @[row:%d, col:%d]  pincer @[row:%d, col:%d]:[row:%d, col:%d] killClue:%d \n",
              rowNumber( pivot ),       colNumber( pivot ),
              rowNumber( firstPincer ), colNumber( firstPincer ),
              rowNumber( otherPincer ), colNumber( otherPincer ), killClueIndex + 1 );
      found = true;
    }
    chainProtection( false, *trio, elementsof( *trio ) );
  }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xyz_WingAnalysis( RowIterator &&rowIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #37 XYZ-Wing Puzzle Solving Technique"
  // 
  // Analysis steps:
  // 1. Find pivot cell of type: XYZ (exact 3 candidates)
  //    The pivot cell has as usual three houses: RH,CH.BH (row, column, box)
  // 2. Find first pincer cell in form: XY (exact 2 candidates) within first (of three) pivot houses
  // 3. Find other pincer cell in form: XZ (exact 2 candidates) within other (of remaining two) pivot houses,
  //    both pincer cells must be in different houses
  // 4. The X-candidate (shared by all 3 XYZ-Wing cells) can be crossed out from the kill zone, which is
  //    defined by every cell sharing a house with the pivot cell and both XYZ-Wing pincers
  //    and are NOT these 3 wing cells itself.
  // 
  //   It makes no sense to iterate row by row _AND_THEN_ column by column when looking for w-Wing constellations.
  bool found = false;
  Cell *trio[3] {};

  for (size_t row = 0; row < SUDOKU_GRID_SIZE; row++, rowIterator.nextHouse()) 
    for (size_t col = 0; col < SUDOKU_GRID_SIZE; col++) {
      const Cell *pivot = trio[0] = &rowIterator[col];
      if (pivot->solved() || pivot->candidateCount() != 3)
        continue;

      BoxIterator firstPincerIterator { pivot->boxOfCellsStart_ };
      for (int box = 0; box < SUDOKU_GRID_SIZE; box++) {
        const Cell *firstPincer = trio[1] = &firstPincerIterator[box];
        if (firstPincer->solved() || firstPincer->candidateCount() != 2 || firstPincer->candidateMisfitsInto( pivot ))
          continue;

        for (auto [firstPincerCandidate, i] = std::pair( -1, 0 ); i < 2; i++) {          // through both candidates
          // take common candidate, add candidates in stages from firstPincer and look for otherPincer
          firstPincerCandidate = firstPincer->candidateNextIndex( firstPincerCandidate );
          auto otherPincerForm = (pivot->candidate_ - firstPincer->candidate_) + firstPincerCandidate;
          if (firstPincer->occupyingDifferentRow( pivot ))
            if (xyz_WingAnalysis( RowIterator { pivot->rowOfCellsStart_ }, &trio, otherPincerForm ))
              found = true;
          if (firstPincer->occupyingDifferentCol( pivot ))
            if (xyz_WingAnalysis( ColIterator { pivot->colOfCellsStart_ }, &trio, otherPincerForm ))
              found = true;
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::xyz_Wing( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (xyz_WingAnalysis( RowIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::skyscraperInvalidConstellation( const CellFishDuo &pair ) {
  return pair.fish_[0]->sharingSameBox( pair.fish_[1] );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::skyscraperInvalidConstellation( const CellFishDuo &head, const CellFishDuo &tail ) {
  // return true, if constellation is invalid for skyscraper:
  // in rowMode:  when corners from BOTH pairs are in the same columns OR
  //              when corners from both pairs are in different columns OR
  //              when pair's corners are NOT in same two columns of boxes 

  if (head.rowMode_) {
    assert( tail.rowMode_ );
    if (head.fish_[0]->colOfCellsStart_ == tail.fish_[0]->colOfCellsStart_ &&
        head.fish_[1]->colOfCellsStart_ == tail.fish_[1]->colOfCellsStart_
        ||
        head.fish_[0]->colOfCellsStart_ != tail.fish_[0]->colOfCellsStart_ &&
        head.fish_[1]->colOfCellsStart_ != tail.fish_[1]->colOfCellsStart_
        ||
        head.fish_[0]->boxOfCellsStart_->colOfCellsStart_ != tail.fish_[0]->boxOfCellsStart_->colOfCellsStart_  ||
        head.fish_[1]->boxOfCellsStart_->colOfCellsStart_ != tail.fish_[1]->boxOfCellsStart_->colOfCellsStart_)
      return true;

    assert( head.fish_[0]->colOfCellsStart_                   == tail.fish_[0]->colOfCellsStart_ &&
            head.fish_[1]->boxOfCellsStart_->colOfCellsStart_ == tail.fish_[1]->boxOfCellsStart_->colOfCellsStart_
            ||
            head.fish_[0]->boxOfCellsStart_->colOfCellsStart_ == tail.fish_[0]->boxOfCellsStart_->colOfCellsStart_ &&
            head.fish_[1]->colOfCellsStart_                   == tail.fish_[1]->colOfCellsStart_ );
  } else {
    assert( tail.rowMode_ == false );
    if (head.fish_[0]->rowOfCellsStart_ == tail.fish_[0]->rowOfCellsStart_ &&
        head.fish_[1]->rowOfCellsStart_ == tail.fish_[1]->rowOfCellsStart_
        ||
        head.fish_[0]->rowOfCellsStart_ != tail.fish_[0]->rowOfCellsStart_ &&
        head.fish_[1]->rowOfCellsStart_ != tail.fish_[1]->rowOfCellsStart_
        ||
        head.fish_[0]->boxOfCellsStart_->rowOfCellsStart_ != tail.fish_[0]->boxOfCellsStart_->rowOfCellsStart_  ||
        head.fish_[1]->boxOfCellsStart_->rowOfCellsStart_ != tail.fish_[1]->boxOfCellsStart_->rowOfCellsStart_)
      return true;

    assert( head.fish_[0]->rowOfCellsStart_                   == tail.fish_[0]->rowOfCellsStart_ &&
            head.fish_[1]->boxOfCellsStart_->rowOfCellsStart_ == tail.fish_[1]->boxOfCellsStart_->rowOfCellsStart_
            ||
            head.fish_[0]->boxOfCellsStart_->rowOfCellsStart_ == tail.fish_[0]->boxOfCellsStart_->rowOfCellsStart_ &&
            head.fish_[1]->rowOfCellsStart_                   == tail.fish_[1]->rowOfCellsStart_ );
  }

  return false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::skyscraperCrossOut( Cell *head, Cell *tail, bool rowMode, int clueIndex ) {
  bool found = false;

  BoxIterator headBoxIterator { head->boxOfCellsStart_ };
  for (int i = 0; i < SUDOKU_BOX_SIZE * SUDOKU_BOX_SIZE; i++) {
    Cell &boxCell = headBoxIterator[i];
    if (rowMode != false && boxCell.colOfCellsStart_ == tail->colOfCellsStart_  ||
        rowMode == false && boxCell.rowOfCellsStart_ == tail->rowOfCellsStart_)          // kill zone?
      if (boxCell.candidatePossible( clueIndex )) {
        printf( "skyscraper cross out clue:%d @[row:%d, col:%d]\n",
                clueIndex + 1, rowNumber( &boxCell ), colNumber( &boxCell ) );
        boxCell.candidateCrossOut( clueIndex );
        found = true;
      }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::skyscraperCrossOut( CellFishDuo *head, CellFishDuo *tail, int clueIndex ) {
  bool found       = false;
  bool rowMode     = head->rowMode_;
  int  cornerIndex = (rowMode != false && head->fish_[0]->colOfCellsStart_ != tail->fish_[0]->colOfCellsStart_
                      ||
                      rowMode == false && head->fish_[0]->rowOfCellsStart_ != tail->fish_[0]->rowOfCellsStart_) ? 0 : 1;
  //if (rowMode)
  //  printf( "skyscraper [rows: %d, %d], [head cols: %d, %d], [tail cols: %d, %d]\n",
  //          rowNumber( head->fish_[0] ), rowNumber( tail->fish_[0] ),
  //          colNumber( head->fish_[0] ), colNumber( head->fish_[1] ),
  //          colNumber( tail->fish_[0] ), colNumber( tail->fish_[1] ) );
  //else
  //  printf( "skyscraper [cols: %d, %d], [head rows: %d, %d], [tail rows: %d, %d]\n",
  //          colNumber( head->fish_[0] ), colNumber( tail->fish_[0] ),
  //          rowNumber( head->fish_[0] ), rowNumber( head->fish_[1] ),
  //          rowNumber( tail->fish_[0] ), rowNumber( tail->fish_[1] ) );

  if (skyscraperCrossOut( head->fish_[cornerIndex], tail->fish_[cornerIndex], rowMode, clueIndex ))
    found = true;
  if (skyscraperCrossOut( tail->fish_[cornerIndex], head->fish_[cornerIndex], rowMode, clueIndex ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::skyscraperAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: intermediate/advanced. Description:
  //   YouTube: "dxSudoku #27 Skyscraper Puzzle Solving Technique"
  // 
  //   Skyscraper constellation requires two pairs of a particular candidate.
  //   Each skyscraper's cell has to be in its own box.
  //   The cells of both pairs have to line up at one pair's end and have to NOT line up at the opposite end.
  //   The second pair cells have to be in boxes, which are below(!) those from the first pair.
  //   The skyscraper's kill zone is determined by the cells that do NOT line up and consists of cells that
  //   share the same house together (the intersection of one cell's box with column/row of the other cell).
  //
  //   Same algorithm can be iterated row by row and then column by column.
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, houseIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      houseIterator[col].mergeCellPossibilities( &mergedSum[row] );

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int headRow = 0; headRow < elementsof( mergedSum ) - SUDOKU_BOX_SIZE; headRow++) {
      if (CellFishDuo::invalidClueCount( mergedSum[headRow].possibilities_[clueIndex].count_ ))
        continue;                                          // skyscraper requires PAIRS of potential clues
      CellFishDuo skyHead { houseIterator, rowMode, mergedSum[headRow], clueIndex };
      if (skyscraperInvalidConstellation( skyHead ))
        continue;
      int tailRow = headRow + SUDOKU_BOX_SIZE - headRow % SUDOKU_BOX_SIZE;     // start row of box below
      for (; tailRow < elementsof(mergedSum); tailRow++) {
        if (CellFishDuo::invalidClueCount( mergedSum[tailRow].possibilities_[clueIndex].count_ ))
          continue;
        CellFishDuo skyTail { houseIterator, rowMode, mergedSum[tailRow], clueIndex };
        if (skyscraperInvalidConstellation( skyTail )  ||  skyscraperInvalidConstellation( skyHead, skyTail ))
          continue;

        if (skyscraperCrossOut( &skyHead, &skyTail, clueIndex ))
          found = true;
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::skyscraper( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (skyscraperAnalysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (skyscraperAnalysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::twoStringKiteCrossOut( const Cell *firstEP, const Cell *otherEP, int killCandidateIndex ) {
  bool found = false;

  if (intersectionOfCells( firstEP, otherEP )->candidateCrossOutInProtectedCell( killCandidateIndex )  ||
      intersectionOfCells( otherEP, firstEP )->candidateCrossOutInProtectedCell( killCandidateIndex ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::twoStringKiteInvalidConstellation( const CellFishDuo &xWingHalf ) {
  return xWingHalf.fish_[0]->sharingSameBox( xWingHalf.fish_[1] );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::twoStringKiteInvalidConstellation( const CellFishDuo &xWingHead, const CellFishDuo &xWingTail ) {

  if (twoStringKiteInvalidConstellation( xWingTail ))
    return true;
  if (xWingHead.fish_[0] == xWingTail.fish_[0]  ||  xWingHead.fish_[0] == xWingTail.fish_[1]  ||
      xWingHead.fish_[1] == xWingTail.fish_[0]  ||  xWingHead.fish_[1] == xWingTail.fish_[1] )
    return true;

  CellFishQuartet quartet {};

  quartet.fish_[0] = xWingHead.fish_[0];
  quartet.fish_[1] = xWingHead.fish_[1];
  quartet.fish_[2] = xWingTail.fish_[0];
  quartet.fish_[3] = xWingTail.fish_[1];

  if (quartet.occupiedBoxesCount() != 3)
    return true;

  return false;              // 2-String Kite occupies 3 boxes
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::twoStringKiteAnalysis( ColIterator &&colIterator, const CellFishDuo &xWingHead, int clueIndex ) {
  bool found = false;

  for (int tailCol = 0; tailCol < SUDOKU_GRID_SIZE; tailCol++, colIterator.nextHouse()) {
    if (candidateCountInHouse( colIterator, clueIndex ) != 2)
      continue;
    CellFishDuo xWingTail { colIterator, false, clueIndex };
    if (twoStringKiteInvalidConstellation( xWingHead, xWingTail ))
      continue;

    int headEndPointIndex = (xWingHead.fish_[0]->occupyingDifferentBox( xWingTail.fish_[0] ) &&
                             xWingHead.fish_[0]->occupyingDifferentBox( xWingTail.fish_[1] )) ? 0 : 1;

    int tailEndPointIndex = (xWingTail.fish_[0]->occupyingDifferentBox( xWingHead.fish_[0] ) &&
                             xWingTail.fish_[0]->occupyingDifferentBox( xWingHead.fish_[1] )) ? 0 : 1;
    
    if (twoStringKiteCrossOut( xWingHead.fish_[headEndPointIndex], xWingTail.fish_[tailEndPointIndex], clueIndex )) {
      printf( "2-String Kite clue:%d @[row:%d, col:%d]:[row:%d, col:%d]:[row:%d, col:%d]:[row:%d, col:%d]\n",
              clueIndex + 1,
              rowNumber( xWingHead.fish_[0] ), colNumber( xWingHead.fish_[0] ),
              rowNumber( xWingHead.fish_[1] ), colNumber( xWingHead.fish_[1] ),
              rowNumber( xWingTail.fish_[0] ), colNumber( xWingTail.fish_[0] ),
              rowNumber( xWingTail.fish_[1] ), colNumber( xWingTail.fish_[1] ) );
      found = true;
    }
  }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::twoStringKiteAnalysis( RowIterator &&rowIterator ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #31 2-String Kite Puzzle Solving Technique"
  //   YouTube: "dxSudoku #74 2-String Kite with Links and Chaining
  //
  // Analysis steps:
  // 1. Search rowwise for a first half of an X-Wing
  // 2. Search columnwise for a second half of an X-Wing (same candidate as in first X-Wing half)
  //    Each X-Wing half must occupy different boxes.
  // 3. 2-String Kite occurs if both halves occupy together exactly 3 boxes.
  // 4. The X-Wing-candidate can be crossed out from a single cell kill zone, which is defined by
  //    the intersection of X-Wing cells, which do NOT share the same box.
  //
  //   It makes no sense to iterate row by row _AND_THEN_ column by column in case of 2-String Kite.
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, rowIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      rowIterator[col].mergeCellPossibilities( &mergedSum[row] );

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int headRow = 0; headRow < elementsof( mergedSum ) - SUDOKU_BOX_SIZE; headRow++) {
      if (CellFishDuo::invalidClueCount( mergedSum[headRow].possibilities_[clueIndex].count_ ))
        continue;                                          // 2-String Kite requires PAIRS of potential clues
      CellFishDuo xWingHead { rowIterator, true, mergedSum[headRow], clueIndex };
      if (twoStringKiteInvalidConstellation( xWingHead ))
        continue;

      if (twoStringKiteAnalysis( ColIterator { xWingHead.fish_[0]->rowOfCellsStart_->colOfCellsStart_ }, xWingHead, clueIndex))
        found = true;
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::twoStringKite( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (twoStringKiteAnalysis( RowIterator { (*grid_).data() } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::emptyRectangleInvalidConstellation( const CellFishDuo &pair ) {
  return (pair.fish_[0]->boxOfCellsStart_ == pair.fish_[1]->boxOfCellsStart_);
}  // -------------------------------------------------------------------------------------------

SudokuGrid::Cell *SudokuGrid::emptyRectangleKillZone( BoxIterator *boxOfCells, const CellFishDuo &pair, int clueIndex ) {
  // return killZone (single cell), if valid AB-switch was found and clueIndex candidate is possible to be crossed out
  // 
  // AB-switch happens when digit's candidates lie exclusively in one row and one column within
  // a box (two perpendicular axes).
  // One of the axes (axisA) has to line up with column (or row) of head.fish_[0] cell.
  // The intersection of the other axis (B) and the column (or row) of pair.fish_[1] cell is the desired kill zone.
  // In rowMode == true the column information of pair cells is relevant else the row information.
  // To get a valid AB-switch at least one candidate cell is needed on axis_A, remaining candidate cells in
  // the boxOfCells have to lie exclusively on ONE axis (axisB).

  const bool &rowMode = pair.rowMode_;
  bool        axisA_found = false;
  Cell       *axisB = nullptr;
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    if ((*boxOfCells)[i].candidatePossible( clueIndex ))
      if (rowMode != false && (*boxOfCells)[i].colOfCellsStart_ == pair.fish_[0]->colOfCellsStart_  ||
          rowMode == false && (*boxOfCells)[i].rowOfCellsStart_ == pair.fish_[0]->rowOfCellsStart_)
        axisA_found = true;
      else if (axisB == nullptr)
        axisB = (rowMode) ? (*boxOfCells)[i].rowOfCellsStart_ : (*boxOfCells)[i].colOfCellsStart_;
      else if (axisB != ((rowMode) ? (*boxOfCells)[i].rowOfCellsStart_ : (*boxOfCells)[i].colOfCellsStart_))
        return nullptr;                          // candidates for axisB occupies different axes

  if (axisA_found == false  ||  axisB == nullptr)
    return nullptr;                              // invalid AB-switch group

  // found valid AB-switch!
  Cell *killCell = (rowMode) ? intersectionOfCells( axisB, pair.fish_[1] )
                             : intersectionOfCells( pair.fish_[1], axisB );
  if (killCell->candidateImpossible( clueIndex ))
    return nullptr;                              // killZone cell candidate is already out
  
  printf( "empty rectangle clue:%d (mode:%s) (AB-switch in box:%d ",
          clueIndex + 1, (rowMode) ? "rows" : "cols", boxNumber( (*boxOfCells).cell_ ) );
  if (rowMode)
    printf( "axisA-col:%d axisB-row:%d", colNumber( pair.fish_[0] ), rowNumber( axisB ) );
  else
    printf( "axisA-row:%d axisB-col:%d", rowNumber( pair.fish_[0] ), colNumber( axisB ) );
  printf( ") cross out @[row:%d, col:%d]\n", rowNumber( killCell ), colNumber( killCell ) );

  return killCell;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::emptyRectangleCrossOut( Cell *(*killZone)[2 * (SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE - 1)], int clueIndex ) {
  bool found = false;

  for (auto killItem : *killZone)
    if (killItem != nullptr  &&  killItem->candidatePossible( clueIndex )) {
      killItem->candidateCrossOut( clueIndex );
      found = true;
      printf( "emptyRectangleCrossOut() clue:%d @[row:%d, col:%d]\n",
              clueIndex + 1, rowNumber( killItem ), colNumber( killItem ) );
    }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::emptyRectangleAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #29 Empty Rectangle Puzzle Solving Technique"
  //   https://www.sudoku.org.pl/pustyprostokat.html
  //
  //   Empty rectangle constellation begins with a pair of a particular candidate located in different boxes.
  //   Then a group called AB-switch (located in a separated box) need to be conjugated with one of the pair cells.
  //   AB-switch happens when digit's candidates lie exclusively in one row and one column (two axes) of a box.
  //   AB-switch group contains between 2 and 5 cells with a particular candidate.
  //   One axis (A) (row or column) points at one of the pair cells, then the other axis (B) points at a potential
  //   kill zone (that is a single cell only).
  //
  //   Same algorithm to find the starting pair of cells can be iterated row by row and then column by column.
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, houseIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      houseIterator[col].mergeCellPossibilities( &mergedSum[row] );

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int pairRow = 0; pairRow < elementsof( mergedSum ); pairRow++) {
      if (CellFishDuo::invalidClueCount( mergedSum[pairRow].possibilities_[clueIndex].count_ ))
        continue;                                          // empty rectangle requires PAIRS of a potential clue
      CellFishDuo pair { houseIterator, rowMode, mergedSum[pairRow], clueIndex };
      if (emptyRectangleInvalidConstellation( pair ))
        continue;                                          // corners are NOT in different boxes

      Cell *killZone[2 * (SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE - 1)] = { nullptr };
      Cell **killZoneIterator = killZone;

      for (int i = 0; i < elementsof( pair.fish_ ); i++) {
        if (i > 0)
          std::rotate( pair.fish_, &pair.fish_[1], &pair.fish_[elementsof( pair.fish_ )]);
        BoxIterator boxOfCells { (rowMode) ? pair.fish_[0]->colOfCellsStart_->boxOfCellsStart_
                                           : pair.fish_[0]->rowOfCellsStart_->boxOfCellsStart_ };
        for (int box = 0; box < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; box++, boxOfCells.nextChuteBox( !rowMode )) {
          if (boxOfCells == pair.fish_[0]->boxOfCellsStart_)
            continue;                                      /* skip own box */
          //if (rowMode)
          //  printf( "empty rectangle: clues:%d pair in row:%d [col:%d : col:%d]\n", clueIndex + 1,
          //          rowNumber( pair.fish_[0] ), colNumber( pair.fish_[0] ), colNumber( pair.fish_[1] ) );
          //else
          //  printf( "empty rectangle: clues:%d pair in col:%d [row:%d : row:%d]\n", clueIndex + 1,
          //          colNumber( pair.fish_[0] ), rowNumber( pair.fish_[0] ), rowNumber( pair.fish_[1] ) );
          //
          //printf( "boxOfCells analysis @%d (box iteration %d) (boxNumer:%d)\n",
          //        (int)(boxOfCells.cell_ - boxOfCells.cell_->rowOfCellsStart_->colOfCellsStart_),
          //        box, boxNumber( boxOfCells.cell_ ) );
          *killZoneIterator++ = emptyRectangleKillZone( &boxOfCells, pair, clueIndex );
        }
      }
      if (emptyRectangleCrossOut( &killZone, clueIndex ))
        found = true;
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::emptyRectangle( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (emptyRectangleAnalysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (emptyRectangleAnalysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------


bool SudokuGrid::finnedX_WingAnalysis( const CellFishDuo &xWingHead, const Cell *third, bool rowMode, int clueIndex ) {
  bool found = false;

  if (third == xWingHead.fish_[0])
    return false;
  if (third->candidateImpossible( clueIndex )  ||  third->sharingSameBox( xWingHead.fish_[0] ))
    return false;

  auto group = groupOfCells( intersectionOfCells( third, xWingHead.fish_[1], rowMode ), rowMode);
  auto &groupIterator = (rowMode) ? (HouseIterator &) RowIterator { &(*group)[0] }
                                  : (HouseIterator &) ColIterator { &(*group)[0] };
  int countInGroup = candidateCountInGroup( groupIterator, clueIndex );
  if (countInGroup == 1)
    return false;                                // --> usual X-Wing or Skyscraper case

  auto &xWingTailIter = (rowMode) ? (HouseIterator &) RowIterator { third->rowOfCellsStart_ }
                                  : (HouseIterator &) ColIterator { third->colOfCellsStart_ };
  if (candidateCountInHouse( xWingTailIter, clueIndex ) - countInGroup != 1)
    return false;                                // --> candidates beside third cell and group of cells/nodes

  candidateGroupProtection( true, groupIterator );
  if (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( BoxIterator { (*group)[0].boxOfCellsStart_ },
                                                               clueIndex, &xWingHead.fish_[1], 1 )) {
    printf( "Finned/Sashimi X-Wing (%swise) clue:%d "
            "@[row:%d, col:%d]:[row:%d, col:%d]:[row:%d, col:%d] groupStart @[row:%d, col:%d]\n",
            (rowMode) ? "row" : "col", clueIndex + 1,
            rowNumber( xWingHead.fish_[0] ), colNumber( xWingHead.fish_[0] ),
            rowNumber( xWingHead.fish_[1] ), colNumber( xWingHead.fish_[1] ),
            rowNumber( third ),              colNumber( third ),
            rowNumber( (*group) ), colNumber( (*group) ));
    found = true;
  }
  candidateGroupProtection( false, groupIterator );
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::finnedX_WingAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: extreme. Description:
  //   YouTube: "dxSudoku #68 Finned X-Wing Puzzle Solving Technique"
  //   YouTube: "dxSudoku #77 Improved Finned X-Wing Search Algorithm"
  //   YouTube: "dxSudoku #76 Sashimi X-Wing"
  // 
  // Finned X-Wing hits partially a subset of Skyscraper and X-Wing cases,
  // therefore it seems Finned X-Wing should be launched after these methods.
  //
  // Analysis steps:
  // 1. Finned or Sahimi X-Wing consists always of 4 boxes.
  //    Find a half of an X-Wing, where both cells occupies 2 different boxes.
  // 2. Find 3rd cell with the X-Wing-Half candidate, which occurs in (separated) third box, but
  //    lines up with any element of found X-Wing half.
  // 3. Analyse fourth box, which lines up with the 3rd cell and the "remaining free" cell of found X-Wing half.
  //    Finned or Sashimi X-Wing happens if the 3rd cell line contains our X-Wing half candidate in the fourth box only.
  // 4. The X-Wing candidate can be crossed out from the kill zone, which is defined by
  //    the intersection of the fourth box X-Wing cells with the "remaining free" cell of found X-Wing half.
  // 
  //   Same algorithm should be iterated row by row and then column by column.
  bool found = false;
  Possibilities mergedSum[SUDOKU_GRID_SIZE] = {0};

  for (int row = 0; row < elementsof( mergedSum ); row++, houseIterator.nextHouse())
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++)
      houseIterator[col].mergeCellPossibilities( &mergedSum[row] );

  for (int clueIndex = 0; clueIndex < elementsof( mergedSum[0].possibilities_ ); clueIndex++)
    for (int headRow = 0; headRow < elementsof( mergedSum ) - SUDOKU_BOX_SIZE; headRow++) {
      if (CellFishDuo::invalidClueCount( mergedSum[headRow].possibilities_[clueIndex].count_ ))
        continue;                                          // skyscraper requires PAIRS of potential clues
      CellFishDuo xWingHead { houseIterator, rowMode, mergedSum[headRow], clueIndex };
      if (xWingHead.fish_[0]->sharingSameBox( xWingHead.fish_[1] ))
        continue;
      for (int i = 0; i < elementsof( xWingHead.fish_ ); i++) {
        if (i > 0)
          std::rotate( xWingHead.fish_, xWingHead.fish_ + 1, xWingHead.fish_ + elementsof( xWingHead.fish_ ) );

        auto &thirdIterator = (rowMode) ? (HouseIterator &) ColIterator { xWingHead.fish_[0]->colOfCellsStart_ }
                                        : (HouseIterator &) RowIterator { xWingHead.fish_[0]->rowOfCellsStart_ };
        for (int k = 0; k < SUDOKU_GRID_SIZE; k++) {
          const Cell *third = &thirdIterator[k];
          if (finnedX_WingAnalysis( xWingHead, third, rowMode, clueIndex ))
            found = true;
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::finnedX_Wing( void ) {
  bool found = false;

  if (finnedX_WingAnalysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (finnedX_WingAnalysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::finnedSwordfish( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: extreme. Description:
  //   YouTube: "Sudoku Swami Finned & Sashimi Swordfish Unveiled! / Tutorial #13"
  //
  // Analysis steps:
  //   Each base set third of Finned Swordfish have to occupy at least 2 boxes.
  //   Particular partly mismathed base set can be accepted as a "fin joker".
  //   Only ONE box of the fish sequence can contain "fin joker", which:
  //    - provides additional cell or cells in base set line (fin case: CFx CxF CFF  [FCx xCF FCF]  [FxC xFC FFC] )
  //    - provides missed cell in base set line on "wrong" position (skyscraper case: .Fx .xF  [F.x x.F]  [xF. Fx.] )
  //    - provides instead of missed cell in base set line two cells (sashimi case: .FF [F.F] [FF.] )
  //
  //   Swordfish candidate cross out can only occur in the same box as the fins.
  //   The kill zone is defined by the intersection of the fin joker box and one or two fish cover sets.
  //   (Cover sets are perpendicular to base sets.)
  //
  //   Finned Swordfish is a variant of X-Chain with group nodes.
  //
  //   Fish base sets searching algorithm can be iterated row by row and then column by column.


  bool found = false;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::finnedFish( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: extreme.
  //
  // Idea description (based on Finned Swordfish - only ONE box can be finned):
  // 1. Try to find Swordfish, Jellyfish, Starfish, Sixfish, Sevenfish, Octopus as ususal
  // 2. When the collection misses only one base set, then:
  //    - rescan the grid for unused lines (=row/col)
  //    - try to apply the "fin joker" if only one box line (group node) spoils the base set
  //    - try to join finned base set with the ordinary base sets.
  // 3. If found the missing base set, cross out only from finned box (in lines perpendicular to the base set).
  //    Cells forming the "fin joker" have to be protected from crossing out.
  // 
  //    "Fin joker" is a group node within a box line (SUDOKU_BOX_SIZE.
  //    "Fin joker" contains further candidate cell/cells in box line even if the candidate is absent at the desired location.
  //    CF. C.F CFF    'C' candidate possible, 'F' = same candidate possible too and spoiling a base set (Fin case)
  //    .F. ..F        'C' candidate missed for base set, but possible within same group   (Skyscraper case)
  //    .SS            'C' candidate missed for base set, but possible twice in same group (Sashimi case)
  //    -----------
  //    C..            'C' candidate possible on a position desired by a particular base set. (No joker required.)
  //    Position of 'C' may vary within a group on each one of three locations, but principle stays same.
  //
  // Is this idea appliable to X-Wing? This is a good question and need to be more investigated.

  bool found = false;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleAnalysisType1( HouseIterator &houseIterator, CellFishQuartet &rectangle ) {
  // 1st, 2nd, 3rd rectangle corner candidates: XY
  // 4th           rectangle corner candidates: XY+[at least one additional candidate]
  // 
  // Kill zone: 4th cell only, but both candidates (X and Y) can be crossed out.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found = false;
  auto &[first, second] = *reinterpret_cast<const Cell *(*)[2]> (rectangle.fish_);

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    const Cell *third = rectangle.fish_[2] = &houseIterator[i];
    if (third == first || third->solved() || third->candidateCount() != 2 || third->candidateDifferent( first ))
      continue;

    Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second, rectangle.rowMode_ );
    if (fourth->solved() || fourth->candidateCount() < 3 || third->candidateMisfitsInto( fourth ))
      continue;
    if (rectangle.occupiedBoxesCount() != 2)
      continue;

    // unique rectangle type 1 found
    if (fourth->candidateCrossOut( third->candidate_ )) {
      if (rectangle.rowMode_)
        printf( "uniqueRectangle (rowwise) [row:%d col:%d,%d; row:%d col:%d] type1 found @[row:%d, col:%d]",
                rowNumber( first ), colNumber( first ), colNumber( second ),
                rowNumber( third ), colNumber( third ),
                rowNumber( fourth ), colNumber( fourth ) );
      else
        printf( "uniqueRectangle (colwise) [row:%d,%d col:%d; row:%d col:%d] type1 found @[row:%d, col:%d]",
                rowNumber( first ), rowNumber( second ), colNumber( first ),
                rowNumber( third ), colNumber( third ),
                rowNumber( fourth ), colNumber( fourth ) );
      printf( " clues crossed out: %c%c%c%c%c%c%c%c%c\n",
              third->candidatePossible( 0 ) ? '1' : '-', third->candidatePossible( 1 ) ? '2' : '-',
              third->candidatePossible( 2 ) ? '3' : '-', third->candidatePossible( 3 ) ? '4' : '-',
              third->candidatePossible( 4 ) ? '5' : '-', third->candidatePossible( 5 ) ? '6' : '-',
              third->candidatePossible( 6 ) ? '7' : '-', third->candidatePossible( 7 ) ? '8' : '-',
              third->candidatePossible( 8 ) ? '9' : '-' );
      found = true;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleAnalysisType2( HouseIterator &houseIterator, CellFishQuartet &rectangle ) {
  // 1st, 2nd rectangle corner candidates: XY
  // 3rd, 4th rectangle corner candidates: XY+Z (i.e. exactly one additional candidate, same for 3rd and 4th)
  // 
  // Kill zone: cells of row/col and box shared by 3rd _AND_ 4th corner.
  //            Z candidate can be crossed out, except the both corners itself.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found   = false;
  bool rowMode = rectangle.rowMode_;
  auto &[first, second] = *reinterpret_cast<const Cell *(*)[2]> (rectangle.fish_);

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    const Cell *third = rectangle.fish_[2] = &houseIterator[i];
    if (third == first || third->solved() || third->candidateCount() != 3 || first->candidateMisfitsInto( third ))
      continue;
    const Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second, rectangle.rowMode_ );
    if (fourth->solved() || fourth->candidateDifferent( third ) || rectangle.occupiedBoxesCount() != 2)
      continue;
    // unique rectangle type 2 found
    chainProtection( true, &rectangle.fish_[2], 2);        // protect: third, fourth
    bool innerFound = false;
    auto killTraits = fourth->candidate_ - first->candidate_;
    HouseIterator &houseIterator = (rowMode) ? (HouseIterator &) RowIterator { fourth->rowOfCellsStart_ }
                                             : (HouseIterator &) ColIterator { fourth->colOfCellsStart_ };
    if (candidateCrossOutInProtectedHouse( houseIterator, killTraits ))
      found = innerFound = true;
    if (boxNumber( third ) == boxNumber( fourth ) &&
        candidateCrossOutInProtectedHouse( BoxIterator { fourth->boxOfCellsStart_ }, killTraits ))
      found = innerFound = true;
    chainProtection( false, &rectangle.fish_[2], 2);       // unprotect: third, fourth

    if (innerFound)
      printf( "uniqueRectangle (%swise) [row:%d,%d col:%d,%d] type2 found, kill zone: %s:%d & box:%d\n",
              (rowMode) ? "row" : "col",
              ((rowMode) ? rowNumber : colNumber)( first ), ((rowMode) ? rowNumber : colNumber)( fourth ),
              ((rowMode) ? colNumber : rowNumber)( first ), ((rowMode) ? colNumber : rowNumber)( fourth ),
              (rowMode) ? "row" : "col", ((rowMode) ? rowNumber : colNumber)( fourth ),
              (boxNumber( third ) == boxNumber( fourth )) ? boxNumber( fourth ) : -1 );
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleAnalysisType3( HouseIterator &houseIterator, CellFishQuartet &rectangle ) {
  // 1st, 2nd rectangle corner candidates: XY
  // 3rd      rectangle corner candidates: XY+[at least one additional candidate]
  // 4th      rectangle corner candidates: XY+[at least one additional candidate, independently from 3rd corner]
  // 
  // The sum of additional candidates of 3rd and 4th corner cell describe virtual cell candidates.
  // If a house shared by 3rd and 4th corner cell contain a naked constellation (starting with the virtual cell),
  // then candidate values can by crossed out as known for nakedDuo, nakedTrio, nakedQuartet, nakedQuintet, ...
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found = false;
  auto &[first, second] = *reinterpret_cast<const Cell *(*)[2]> (rectangle.fish_);

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *third = rectangle.fish_[2] = &houseIterator[i];
    if (third == first || third->solved() || third->candidateCount() < 3 || first->candidateMisfitsInto( third ))
      continue;
    Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second, rectangle.rowMode_ );
    if (fourth->solved() || fourth->candidateCount() < 3 || first->candidateMisfitsInto( fourth ) ||
        rectangle.occupiedBoxesCount() != 2)
      continue;
    const auto virtualCellTraits = (third->candidate_ - first->candidate_) + (fourth->candidate_ - first->candidate_);
    if (virtualCellTraits.candidateCount() < 2)
      continue;

    if (NakedMethodDuo     {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;
    if (NakedMethodTrio    {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;
    if (NakedMethodQuartet {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;
    if (NakedMethodQuintet {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;
    if (NakedMethodSextet  {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;
    if (NakedMethodSeptet  {}.analysisURtype3( third, fourth, virtualCellTraits ))
      found = true;

    if (found) {
      printf( "uniqueRectangle (%swise) [row:%d,%d col:%d,%d] TYPE3 found, rectBase @[row:%d, col:%d] @[row:%d, col:%d]\n",
              (rectangle.rowMode_) ? "row" : "col",
              rowNumber( first ), rowNumber( fourth ), colNumber( first ), colNumber( fourth ),
              rowNumber( third ), colNumber( third ),  rowNumber( fourth ), colNumber( fourth ) );
      break;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleAnalysisType4( HouseIterator &houseIterator, CellFishQuartet &rectangle ) {
  // 1st, 2nd rectangle corner candidates: XY
  // 3rd      rectangle corner candidates: XY+[at least one additional candidate]
  // 4th      rectangle corner candidates: XY+[at least one additional candidate, independently from 3rd corner]
  // 
  // Kill zone: if X candidate occurs precisely 2 times in row/column shared by 3rd and 4th corner,
  //            then the Y candidate can be crossed out from these both corner cells.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found = false;
  auto &[first, second] = *reinterpret_cast<const Cell *(*)[2]> (rectangle.fish_);

  for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *third = rectangle.fish_[2] = &houseIterator[i];
    if (third == first || third->solved() || third->candidateCount() < 3 || first->candidateMisfitsInto( third ))
      continue;
    Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second, rectangle.rowMode_ );
    if (fourth->solved() || fourth->candidateCount() < 3 || first->candidateMisfitsInto( fourth ) ||
        rectangle.occupiedBoxesCount() != 2)
      continue;

    auto &rectBaseIterator = (rectangle.rowMode_) ? (HouseIterator &) RowIterator { fourth->rowOfCellsStart_ }
                                                  : (HouseIterator &) ColIterator { fourth->colOfCellsStart_ };
    for (auto [k, firstCornerCandidate] = std::pair { 0, -1 }; k < 2; k++) {   // through both XY candidates
      firstCornerCandidate = first->candidateNextIndex( firstCornerCandidate );
      if (candidateCountInHouse( rectBaseIterator, firstCornerCandidate ) > 2)
        continue;

      // unique rectangle type 4 found, because an X pair of candidates
      // occurs precisely twice in row/column shared by third and fourth corner.
      int killClueIndex = first->candidate_ - firstCornerCandidate;            // subtraction produces here only one candidate
      printf( "uniqueRectangle (%swise) [row:%d,%d col:%d,%d] TYPE4 found, clue:%d cross out @[row:%d, col:%d] @[row:%d, col:%d]\n",
              (rectangle.rowMode_) ? "row" : "col",
              rowNumber( first ), rowNumber( fourth ), colNumber( first ), colNumber( fourth ), killClueIndex + 1,
              rowNumber( third ), colNumber( third ),  rowNumber( fourth ), colNumber( fourth ) );

      third->candidateCrossOut( killClueIndex );
      fourth->candidateCrossOut( killClueIndex );
      found = true;
      break;       // third and fourth corners modified!!
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleAnalysis( HouseIterator &&houseIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #45 Unique Rectangle Type 1"
  //   YouTube: "dxSudoku #46 Unique Rectangle Type 2"
  //   YouTube: "dxSudoku #47 Unique Rectangle Type 3"
  //   YouTube: "dxSudoku #48 Unique Rectangle Type 4"
  // 
  //   It makes no sense to iterate row by row _AND_THEN_ column by column in case of type 1.
  //   Same algorithm should be iterated row by row and then column by column in case of type 2, 3, 4.
  bool found = false;
  CellFishQuartet rectangle;
  rectangle.rowMode_ = rowMode;

  for (int row = 0; row < SUDOKU_GRID_SIZE; row++, houseIterator.nextHouse())
    for (int frontCol = 0; frontCol < SUDOKU_GRID_SIZE - 1; frontCol++) {
      const Cell *first = rectangle.fish_[0] = &houseIterator[frontCol];
      if (first->solved() || first->candidateCount() != 2)
        continue;
      for (int rearCol = frontCol + 1; rearCol < SUDOKU_GRID_SIZE; rearCol++) {
        const Cell *second = rectangle.fish_[1] = &houseIterator[rearCol];
        if (second->solved() || second->candidateCount() != 2 || second->candidateDifferent( first ))
          continue;

        for (int i = 0; i < 2; i++) {
          if (i > 0)
            std::rotate( rectangle.fish_, &rectangle.fish_[1], &rectangle.fish_[2] );

          if (rowMode)
            if (uniqueRectangleAnalysisType1( ColIterator { first->colOfCellsStart_ }, rectangle ))
              found = true;

          if (i == 0) {
            HouseIterator &thirdCornerIterator = (rowMode) ? (HouseIterator &) ColIterator { first->colOfCellsStart_ }
                                                           : (HouseIterator &) RowIterator { first->rowOfCellsStart_ };
            if (uniqueRectangleAnalysisType2( thirdCornerIterator, rectangle ))
              found = true;
            if (uniqueRectangleAnalysisType3( thirdCornerIterator, rectangle ))
              found = true;
            if (uniqueRectangleAnalysisType4( thirdCornerIterator, rectangle ))
              found = true;
          }
        }
      }
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleType5CrossOut( const Cell *rectangleCorner, CellFishQuartet &rectangle, int clueIndex ) {
  bool   found = false;
  size_t neighboursQuantity = (rectangle.fish_[3]->candidateCount() < 3) ? 2 : 3;
#if 1
  if (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( BoxIterator { rectangleCorner->boxOfCellsStart_ },
       clueIndex, &rectangle.fish_[1], neighboursQuantity ))
    found = true;
#else
  BoxIterator killBoxOfCells { rectangleCorner->boxOfCellsStart_ };

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &killBoxOfCells[i];
    if (cell->inAllNeighboursSight( const_cast<const Cell **> (rectangle.fish_ + 1), neighboursQuantity ) &&
        cell->candidateCrossOutInProtectedCell( clueIndex )) {
      printf( " @[row:%d, col:%d]", rowNumber( cell ), colNumber( cell ) );
      found = true;
    }
  }
#endif
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleScrutinyType5( RowIterator &rowIterator, CellFishQuartet &rectangle ) {
  // 1st      rectangle          corner candidates: XY
  // 2nd, 3rd rectangle diagonal corner candidates: XYZ (i.e. _ONE_ additional candidate)
  // 4th      rectangle          corner candidates: XY or XYZ (with same additional candidate)
  // 
  // Kill zone: cells of row/col/box shared by (2 or 3) rectangle corners with additional Z candidate.
  //            Z candidate can be crossed out, except from the corners itself.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found = false;

  const Cell *first = rectangle.fish_[0];
  for (int col = 0; col < SUDOKU_GRID_SIZE; col++) {
    const Cell *second = rectangle.fish_[1] = &rowIterator[col];
    if (second == first || second->solved() || second->candidateCount() != 3 || first->candidateMisfitsInto( second ))
      continue;
    ColIterator colOfCells { first->colOfCellsStart_ };
    for (int row = 0; row < SUDOKU_GRID_SIZE; row++) {
      const Cell *third = rectangle.fish_[2] = &colOfCells[row];
      if (third == first || third->solved() || third->candidateCount() != 3 || third->candidateDifferent( second ))
        continue;
      const Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second );
      if (fourth->solved() || fourth->candidateDifferent( first ) && fourth->candidateDifferent( third ))
        continue;
      if (rectangle.occupiedBoxesCount() != 2)
        continue;

      // unique rectangle type 5 found
      int killClueIndex = third->candidate_ - first->candidate_;     // subtraction produces here only one candidate
      rectangle.chainProtection( true );
      assert( boxNumber ( second ) != boxNumber( third ) );
      if (uniqueRectangleType5CrossOut( second, rectangle, killClueIndex ))
        found = true;
      if (uniqueRectangleType5CrossOut( third,  rectangle, killClueIndex ))
        found = true;
      rectangle.chainProtection( false );
      if (found)
        printf( "uniqueRectangle [row:%d,%d col:%d,%d] TYPE5 found, clue:%d crossed out\n",
                rowNumber( first ), rowNumber( fourth ),
                colNumber( first ), colNumber( fourth ), killClueIndex + 1 );
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleScrutinyType6( RowIterator &rowIterator, CellFishQuartet &rectangle ) {
  // 1st, 4th rectangle diagonal corner candidates: XY
  // 2nd, 3rd rectangle diagonal corner candidates: XY+[at least one additional candidate (but both corners _together_)]
  // 
  // Kill zone: 2nd and 3rd diagonal corner cells (diagonal without additional candidate(s) is not relevant).
  //            X candidate can be crossed out, if X candidate occurs precisely 2 times
  //            in row and column located by 2nd and 3rd corner cell.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  bool found = false;

  const Cell *first = rectangle.fish_[0];        // XY candidate
  for (int col = 0; col < SUDOKU_GRID_SIZE  &&  found == false; col++) {
    Cell *second = rectangle.fish_[1] = &rowIterator[col];
    if (second == first || second->solved() || second->candidateCount() < 2 || first->candidateMisfitsInto( second ))
      continue;                                  // not XY+[optionally any additional candidates]
    ColIterator colOfCells { first->colOfCellsStart_ };
    for (int row = 0; row < SUDOKU_GRID_SIZE  &&  found == false; row++) {
      Cell *third = rectangle.fish_[2] = &colOfCells[row];
      if (third == first || third->solved() || third->candidateCount() < 2 || first->candidateMisfitsInto( third ))
        continue;                                // not XY+[optionally any additional candidates]
      if (second->candidateCount() + third->candidateCount() < 5)
        continue;                                // second and third cannot be both XY only
      const Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second );
      if (fourth->solved() || fourth->candidateDifferent( first ))
        continue;                                // different from XY
      if (rectangle.occupiedBoxesCount() != 2)
        continue;

      for (auto [i, killClueIndex] = std::pair {0U, -1}; i < 2; i++) {
        // analyse both candidates from first corner cell
        killClueIndex = first->candidateNextIndex( killClueIndex );
        if (candidateCountInHouse( RowIterator( second->rowOfCellsStart_ ), killClueIndex ) != 2  ||
            candidateCountInHouse( ColIterator( second->colOfCellsStart_ ), killClueIndex ) != 2  ||
            candidateCountInHouse( RowIterator(  third->rowOfCellsStart_ ), killClueIndex ) != 2  ||
            candidateCountInHouse( ColIterator(  third->colOfCellsStart_ ), killClueIndex ) != 2 )
          continue;

        // unique rectangle type 6 found
        printf( "uniqueRectangle [row:%d,%d col:%d,%d] TYPE6 found, clue:%d crossed out:",
                rowNumber( first ), rowNumber( fourth ),colNumber( first ), colNumber( fourth ), killClueIndex + 1 );
        if (second->candidatePossible( killClueIndex )) {
          second->candidateCrossOut( killClueIndex );
          found = true;
          printf( " @[row:%d, col:%d]", rowNumber( second ), colNumber( second ) );
        }
        if (third->candidatePossible( killClueIndex )) {
          third->candidateCrossOut( killClueIndex );
          found = true;
          printf( " @[row:%d, col:%d]", rowNumber( third ), colNumber( third ) );
        }
        printf( "\n" );
        break;     // second and/or third corner modified!
      }
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleScrutinyType7( RowIterator &rowIterator, CellFishQuartet &rectangle ) {
    // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #56 Hidden Rectangle" aka Unique Rectangle Type 7
  //
  // Analysis steps:
  // 1. Find 1st rectangle corner cell:          XY (exact 2 candidates)
  // 2. Find 2nd corner (in same row):           XY+[at least one additional candidate]
  // 3. Find 3rd corner (in same column as 1st): XY+[at least one additional candidate (independently from 2nd corner)]
  // 4. Find 4th corner (junction of 2nd + 3rd): XY+[none or any additional candidate(s), independently from 2nd/3rd]
  // 5. Hidden Rectangle occurs if X candidate happens exactly 2 times in row and column of the 4th corner 
  // 6. Y candidate can be crossed out from the kill zone, which is the 4th corner cell itself.
  // 
  // All rectangle corners must occur over exactly 2 rows, 2 columns and 2 boxes.
  // 
  // It makes no sense to iterate row by row _AND_THEN_ column by column.
  bool found = false;

  const Cell *first = rectangle.fish_[0];        // XY candidate
  for (int col = 0; col < SUDOKU_GRID_SIZE; col++) {
    const Cell *second = rectangle.fish_[1] = &rowIterator[col];
    if (second == first || second->solved() || second->candidateCount() < 3 || first->candidateMisfitsInto( second ))
      continue;                                  // not XY+[at least one additional candidate]
    ColIterator colOfCells { first->colOfCellsStart_ };
    for (int row = 0; row < SUDOKU_GRID_SIZE; row++) {
      const Cell *third = rectangle.fish_[2] = &colOfCells[row];
      if (third == first || third->solved() || third->candidateCount() < 3 || first->candidateMisfitsInto( third ))
        continue;                                // not XY+[at least one additional candidate]
      Cell *fourth = rectangle.fish_[3] = intersectionOfCells( third, second );
      if (fourth->solved() || fourth->candidateCount() < 2 || first->candidateMisfitsInto( fourth ))
        continue;                                // not XY+[none or any additional candidate(s)]
      if (rectangle.occupiedBoxesCount() != 2)
        continue;
      // unique rectangle type 7 found
      for (auto [i, firstCornerCandidate] = std::pair { 0, -1 }; i < 2; i++) {  // through both candidates
        firstCornerCandidate = first->candidateNextIndex( firstCornerCandidate );
        if (candidateCountInHouse( RowIterator { fourth->rowOfCellsStart_ }, firstCornerCandidate ) != 2 ||
            candidateCountInHouse( ColIterator { fourth->colOfCellsStart_ }, firstCornerCandidate ) != 2)
          continue;

        fourth->candidateCrossOut( (first->candidate_ - firstCornerCandidate).candidateFirstIndex() );
        found = true;
        printf( "hiddenRectangle (unique rectangle type7) @[row:%d, col:%d]:[...]:[...]:[row:%d, col:%d] kill clue:%d\n",
                rowNumber( first ), colNumber( first ), rowNumber( fourth ), colNumber( fourth ),
                (first->candidate_ - firstCornerCandidate).candidateFirstIndex() + 1);
        break;     // fourth corner cell modified!
      }
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangleScrutiny( RowIterator &&rowIterator, bool rowMode ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #52 Unique Rectangle Type 5"
  //   YouTube: "dxSudoku #53 Unique Rectangle Type 6"
  //   YouTube: "dxSudoku #56 Hidden Rectangle" aka "Unique Rectangle Type 7"
  // 
  //   It makes no sense to iterate row by row _AND_THEN_ again column by column in case of types 5, 6, 7.
  bool found = false;
  CellFishQuartet rectangle;
  rectangle.rowMode_ = rowMode;

  for (int row = 0; row < SUDOKU_GRID_SIZE; row++, rowIterator.nextHouse())
    for (int frontCol = 0; frontCol < SUDOKU_GRID_SIZE - 1; frontCol++) {
      const Cell *first = rectangle.fish_[0] = &rowIterator[frontCol];
      if (first->solved() || first->candidateCount() != 2)
        continue;
      if (uniqueRectangleScrutinyType5( rowIterator, rectangle ) ||
          uniqueRectangleScrutinyType6( rowIterator, rectangle ) ||
          uniqueRectangleScrutinyType7( rowIterator, rectangle ))
        found = true;
    }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::uniqueRectangle( void ) {
  // return false, if neither clue possibility can be crossed out
  bool found = false;

  if (uniqueRectangleAnalysis( RowIterator { (*grid_).data() }, true ))
    found = true;
  if (uniqueRectangleAnalysis( ColIterator { (*grid_).data() }, false ))
    found = true;

  if (uniqueRectangleScrutiny( RowIterator { (*grid_).data() }, true ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

#if 0
bool SudokuGrid::hiddenRectangle( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #56 Hidden Rectangle" (aka Unique Rectangle Type 7)
  //
  //   It makes no sense to iterate row by row _AND_THEN_ column by column.
  bool found = false;

  //if (hiddenRectangleAnalysis( RowIterator { (*grid_).data() }, true ))
  //  found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::turbotFish( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #36 Turbot Fish Puzzle Solving Technique"
  //            A Turbot Fish is a X-chain of four cells and three links. First and last link must be a Strong Link.
  //            The kill zone is defined by any cell sharing a house with the endpoints in the chaining sequence.
  // 
  //            The Turbot fish is an X-Chain constellation with exactly 4 cells and 3 links. 
  bool found = false;

  return found;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::X_Chain::Link::chainKillZone( X_Chain *chain, const Cell *endpoint, const Cell *strongLinkStart,
                                               HouseIterator &&houseIterator ) {
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &houseIterator[i];
    if (cell->candidateImpossible( clueIndex ))
      continue;
    if (cell == endpoint || cell == strongLinkStart || cell->occupyingDifferentHouse( strongLinkStart ))
      continue;

    assert( cell->unsolved() );
    // intentionally ignore cell->candidateProtected() status, because it's used by x_chain() for chain trackking purpose
    cell->candidateCrossOut( clueIndex );

    // todo! w przysz³oœci: zanotowaæ do wykreœlenia w strukturze *chain, zamiast natychmiastowego wkreœlania
    //       idea: byæ mo¿e w ten sposób uda siê z³o¿yæ wiêcej ³añcuchów x-chain co powinno zwiêkszyæ szanse
    //       na wiêcej eliminacji.

    printf( "chainKillZone( clue:%d ) crossed out @[row:%d, col:%d]\n", clueIndex + 1, rowNumber( cell ), colNumber( cell ) );
    found = true;
  }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::chainKillZone( X_Chain *chain, const Cell *endpoint, const Cell *strongLinkStart ) {
  // look for intersections between all three houses of 'endpoint' and EVERY strong link starting point in links chain
  bool found = false;

  if (chainKillZone( chain, endpoint, strongLinkStart, RowIterator { endpoint->rowOfCellsStart_ } ))
    found = true;
  if (chainKillZone( chain, endpoint, strongLinkStart, ColIterator { endpoint->colOfCellsStart_ } ))
    found = true;
  if (chainKillZone( chain, endpoint, strongLinkStart, BoxIterator { endpoint->boxOfCellsStart_ } ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::chainKillZone( X_Chain *chain, const Cell *endpoint ) {
  // apply kill zone between 'endpoint' cell and EVERY strong link starting point in links chain.
  // current link object must be weak link ending in links chain (and starting point
  // for a strong link leading to 'endpoint' cell)
  bool found = false;
  const int clueIndex = chain->getClueIndex();
  const Link *link = this;
  
  do {
    assert( link->prev_ != nullptr );
    link = link->prev_->prev_;
    printf( "chainKillZone( clue:%d ) [%d,%d :: %d,%d] \n",
            clueIndex + 1, rowNumber( endpoint ), colNumber( endpoint ), rowNumber( link->cell_ ), colNumber( link->cell_ ) );
    if (chainKillZone( chain, endpoint, link->cell_ ))
      found = true;
  } while (link->prev_ != nullptr);

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::chainClosedKillZone( X_Chain *chain, const Cell *linkEnd, const Cell *linkStart ) {
  // kill zone is a house common for weak link endpoints (strong link consists of 2 cells with particular candidate only).
  // x-chain candidate can be crossed out in that house except of weak link ending and starting cells.
  // return true, if any candidate has been crossed out.
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  assert( linkEnd != linkStart );
  assert( linkEnd->sharingSameRow( linkStart ) ||
          linkEnd->sharingSameCol( linkStart ) ||
          linkEnd->sharingSameBox( linkStart ) );

  HouseIterator &houseIt = linkEnd->sharingSameRow( linkStart ) ? (HouseIterator &) RowIterator { linkEnd->rowOfCellsStart_ } :
                           linkEnd->sharingSameCol( linkStart ) ? (HouseIterator &) ColIterator { linkEnd->colOfCellsStart_ } :
                                                                  (HouseIterator &) BoxIterator { linkEnd->boxOfCellsStart_ };

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &houseIt[i];
    if (cell->candidatePossible( clueIndex ) && cell != linkEnd && cell != linkStart) {
      printf( "chainClosedKillZone(): weak link house used as kill zone - clue%d @[row:%d, col:%d] %s crossed out\n",
              clueIndex + 1, rowNumber( cell ), colNumber( cell ),
              cell->candidateProtected() ? "(cell IN CHAIN !!)" : "" );
      if (cell->candidateProtected())
        printf( "" );
      // intentionally ignore cell->candidateProtected() status, because it's used by x_chain() for chain trackking purpose
      cell->candidateCrossOut( clueIndex );
      found = true;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::chainClosedKillZone( X_Chain *chain /* const Cell *endpoint */) {
  //  Closed x-chain end (last strong link end (= 'endpoint') is in same house as chain root (= chain beginning).
  //  It means: every chained weak link becomes effectively a strong link (Either-Or link), i.e. chain candidates
  //  can be removed from houses defined by all weak links (except of the weak link cells itself).
  //  The loop closing weak link is NOT eliminated here, because it is a job of previously called chainKillZone() method.
  //  Kill zone mechanism explanation: Sudoku Swami - Tutorial #27 mark 35:25, mark 36:50; tutorial #42 mark 20:20
  //
  //  Current link object must be weak link ending in chain (and starting point for a strong link leading to 'endpoint'
  //  cell, which shares a house with chain root).
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  // 'endpoint' seems to be obsolete!

  for (const Link *link = this; link->prev_ != nullptr; link = link->prev_->prev_) {
    assert( link->prev_->prev_ != nullptr );
    printf( "chainClosedKillZone( clue:%d ) link [%d,%d <-- %d,%d] used as kill zone\n",
            clueIndex + 1,
            rowNumber( link->cell_ ),        colNumber( link->cell_ ),
            rowNumber( link->prev_->cell_ ), colNumber( link->prev_->cell_ ) );
    if (chainClosedKillZone( chain, link->cell_, link->prev_->cell_ ))
      found = true;
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::isStrongLinkStartInChain( const Cell *cell ) {
  // return true, if 'cell' is a starting point of already CHAINED strong link.
  // Current link object is assumed to be a chained weak link ending and a potential starting point
  // for a strong link leading to 'cell' (i.e.: ... --> 'this->cell_' ==> 'cell' )  .
  // 'cell' must belong to links in chain.

  assert( cell != this->cell_ );                 // check assumptions
  assert( this->cell_->candidateProtected() );
  assert( cell->candidateProtected() );          // cell must belong to links in chain

  for (const Link *link = this->prev_; link != nullptr; ) {
    if (link->cell_ == cell)
      return false;                              // found cell_ is a strong link ENDING side

    assert( link->prev_ != nullptr );
    link = link->prev_;

    if (link->cell_ == cell)                     // found cell_ is a strong link STARTING side
      return true;

    link = link->prev_;
  }
  assert( false );                               // 'cell' not in chain
  return false;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_Chain::Link::printChain( const char *title, const Cell *cell, const Link *link, const char *msg ) {

  assert( cell != nullptr && link != nullptr && link->prev_ != nullptr );
  printf( "%s", title);
  printf( "chain in reversed order: [%d,%d] <== [%d,%d]",
          rowNumber( cell ), colNumber( cell ), rowNumber( link->cell_ ), colNumber( link->cell_ ) );

  for (link = link->prev_; link != nullptr ; link = link->prev_->prev_) {
    printf( " <-- [%d,%d] <== [%d,%d]", 
            rowNumber( link->cell_ ),        colNumber( link->cell_ ),
            rowNumber( link->prev_->cell_ ), colNumber( link->prev_->cell_ ) );
  }
  printf( "%s", msg);

}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::addStrongLink( X_Chain *chain, HouseIterator &&houseIterator ) {
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  if (candidateCountInHouse( houseIterator, clueIndex ) != 2)
    return false;                      // strong link unavailable

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &houseIterator[i];
    if (cell->candidateImpossible( clueIndex ) || cell == cell_)
      continue;                        // ignore invalid cells and strong link start (== chain ending stored in cell_)

    printf( "addStrongLink(): weak/strong link @[%d,%d --> %d,%d ==> %d,%d]\n",
            rowNumber( prev_->cell_ ), colNumber( prev_->cell_ ),
            rowNumber( cell_ ),        colNumber( cell_ ),
            rowNumber( cell ),         colNumber( cell ) );
    if (cell == chain->getRoot()) {    // special case: strong link end hits chain's start ('nice loop')
      printChain( "addStrongLink(): chain closed by strong link - root cell can be established\n", cell, this, "\n" );
      establishClueAndCrossItOutInRowColBox( cell, clueIndex );
      found = true;
    } else if (cell->candidateProtected()) {     // cell already in chain?
      printChain( "-------------------- but endpoint already in chain and protected.\n", cell, this, "\n" );
      if (isStrongLinkStartInChain( cell )) {    // 'nice loop' as a subsequence of found x-chain
        printf( "addStrongLink(): closed loop ==> cell [%d,%d] SOLVED (%d) !!\n",
                rowNumber( cell ), colNumber( cell ), clueIndex + 1 );
        establishClueAndCrossItOutInRowColBox( cell, clueIndex );
        found = true;
      } else if (chainKillZone( chain, cell ))
        found = true;
    } else {
      //printf( "---- enhanced kill zone mechanism:\n" );
      if (chainKillZone( chain, cell ))            // enhanced kill zone mechanism
        found = true;

      if (cell->sharingSameHouse( chain->getRoot() )) {
        printChain( "addStrongLink(): ---- chain closed by weak link.\n", cell, this,
                    "\nall weak links in loop are futher potential kill zones\n" );
        if (chainClosedKillZone( chain /*, cell */ ))
          found = true;
      }
      assert( cell->candidateUnprotected() );
      cell->candidateProtection( true );
      Link link { cell, this };                    // append new link to the chain
      if (link.addWeakLink( chain ))
        found = true;
      cell->candidateProtection( false );
    }
    break;
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::addStrongLink( X_Chain *chain ) {
  // append any possible "strong link", apply kill zone and go recursively further 
  bool found = false;

  // try to find each possible strong link outgoing from "this->cell_"
  // if cross out successful in kill zone, found = true;
  // add found "strong link" cell 
  // try to find next "weak link" cell

  if (cell_->occupyingDifferentRow( prev_->cell_ ))
    if (addStrongLink( chain, RowIterator { cell_->rowOfCellsStart_ } ))
      found = true;
  if (cell_->occupyingDifferentCol( prev_->cell_ ))
    if (addStrongLink( chain, ColIterator { cell_->colOfCellsStart_ } ))
      found = true;
  if (cell_->occupyingDifferentBox( prev_->cell_ ))
    if (addStrongLink( chain, BoxIterator { cell_->boxOfCellsStart_ } ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::addWeakLink( X_Chain *chain, HouseIterator &&houseIterator ) {
  // append any possible "weak link", try with all potential cells in range of given houseIterator
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  // loop: add each possible "weak link" cell {
  //   znaj¹c kierunek ostatniego strong link in chain (prev_->cell_ ---> cell_)
  //   wyznacz kolejno mo¿liwe "weak links" wychodz¹ce z cell_ (i nietrafiaj¹ce w komórki ju¿ u¿yte w ³añcuchu)
  // 
  //   ka¿d¹ z takich œwie¿ych komórek (weakLinkCell) umieœæ w œwie¿ym ogniwie i dodaj do ³añcucha:
  //   X_Chain link { weakLinkCell, this };
  //   if (link.addStrongLink( clueIndex ))
  //     found = true;
  // }

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &houseIterator[i];
    if (cell->candidatePossible( clueIndex ) && cell->candidateUnprotected()) {
      cell->candidateProtection( true );
      Link link { cell, this };                  // append new link to the chain
      //printf( "weakLink --> [%d,%d]\n", rowNumber( cell ), colNumber( cell ) );
      if (link.addStrongLink( chain ))
        found = true;
      cell->candidateProtection( false );
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::Link::addWeakLink( X_Chain *chain ) {
  // append any possible "weak link", skip searching in the house occupied by last strong link
  bool found = false;

  //printf( "addWeakLink( clue:%d ) Either-Or link: @[%d,%d ==> %d,%d]\n", chain->getClueIndex() + 1,
  //        rowNumber( prev_->cell_ ), colNumber( prev_->cell_ ),
  //        rowNumber( cell_ ),        colNumber( cell_ ) );

  if (cell_->occupyingDifferentRow( prev_->cell_ ))
    if (addWeakLink( chain, RowIterator { cell_->rowOfCellsStart_ } ))
      found = true;
  if (cell_->occupyingDifferentCol( prev_->cell_ ))
    if (addWeakLink( chain, ColIterator { cell_->colOfCellsStart_ } ))
      found = true;
  if (cell_->occupyingDifferentBox( prev_->cell_ ))
    if (addWeakLink( chain, BoxIterator { cell_->boxOfCellsStart_ } ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::seedStrongLink( HouseIterator &houseIterator ) {
  size_t candidateCount = 0;
  const int clueIndex = getClueIndex();

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (houseIterator[i].candidatePossible( clueIndex ))
      if (candidateCount < elementsof( seedLink ))
        seedLink[candidateCount++].cell_ = &houseIterator[i];
      else
        candidateCount += 1;

  return (candidateCount == elementsof( seedLink )) ? true : false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_Chain::seedChain( HouseIterator &&houseIterator, int clueIndex ) {
  bool found = false;

  X_Chain chain { clueIndex };

  for (size_t house = 0; house < SUDOKU_GRID_SIZE; house++, houseIterator.nextHouse()) {
    if (chain.seedStrongLink( houseIterator )) {
    #if 0
      printf( "startChain( clue:%d )  Either-Or link: @[%d,%d --> %d,%d]\n", clueIndex + 1,
              rowNumber( chain.seedLink[0].cell_ ),        colNumber( chain.seedLink[0].cell_ ),
              rowNumber( chain.seedLink[0].prev_->cell_ ), colNumber( chain.seedLink[0].prev_->cell_ ));
    #endif
      assert( chain.seedLink[0].prev_->prev_ == nullptr  &&  chain.seedLink[1].prev_ == nullptr );
      chain.seedLink[0].cell_->candidateProtection( true ), chain.seedLink[1].cell_->candidateProtection( true );
      //printf( "F:startLink( clue:%d ) Either-Or link: @[%d,%d ==> %d,%d]\n", clueIndex + 1,
      //        rowNumber( chain.seedLink[0].prev_->cell_ ), colNumber( chain.seedLink[0].prev_->cell_ ),
      //        rowNumber( chain.seedLink[0].cell_ ),        colNumber( chain.seedLink[0].cell_ ) );
      if (chain.seedLink[0].addWeakLink( &chain ))         // strong link head as chain start
        found = true;

      std::swap( chain.seedLink[0].cell_, chain.seedLink[1].cell_ );
      assert( chain.seedLink[0].prev_->prev_ == nullptr  &&  chain.seedLink[1].prev_ == nullptr );
      //printf( "B:startLink( clue:%d ) Either-Or link: @[%d,%d ==> %d,%d]\n", clueIndex + 1,
      //        rowNumber( chain.seedLink[0].prev_->cell_ ), colNumber( chain.seedLink[0].prev_->cell_ ),
      //        rowNumber( chain.seedLink[0].cell_ ),        colNumber( chain.seedLink[0].cell_ ) );
      if (chain.seedLink[0].addWeakLink( &chain ))         // strong link tail as chain start
        found = true;
      chain.seedLink[0].cell_->candidateProtection( false ), chain.seedLink[1].cell_->candidateProtection( false );
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::x_Chain( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #41 X-Chain Puzzle Solving Technique"
  //   YouTube: "dxSudoku #73 Improved X-Chain Search Algorithm"
  //   YouTube: "dxSudoku #75 Empty Rectangle with Group Nodes and Group Links (marks: 4:30, 5:34)"
  //   YouTube: "dxSudoku #77 Improved Finned X-Wing Search Algorithm" (mark 1:30)
  //   YouTube: "dxSudoku #100 X-Nodes with X-Chains"
  //   YouTube: "dxSudoku #94 X-Chains with Nice Loops"
  //   YouTube: "Sudoku Swami #27 (marks: 19:20, 22:15, 23:08, 23:48, 
  //                               24:12, 24:34 (= cannibalistic PLUS one elimination NOT between chain ends)"
  //
  //   Variants of X-Chain: X-Wing, Skyscraper, 2-String Kite, Empty Rectangle, Turbot Fish, Finned X-Wing, Sashimi X-Wing.
  //
  // Analysis steps:
  // 1. Find two cells with X candidate building a strong link: only 2 cells with particular candidate in a house (row/col/box).
  // 2. Call one end of this strong link as X-chain root and the opposite end as non-root.
  // 3. To build an X-chain a PAIR of further X-candidate cells is needed.
  //    First required pair element can be any (X candidate) cell in sight of the non-root cell (except the root cell itself).
  //    The connection between non-root cell and first pair element forms a "weak link".
  //    Second required pair element must build a "strong link" starting from the first pair element.
  //    All possible "weak link / strong link" pairs should be checked.
  // 4. If described PAIR is found, append it to the X-chain, apply the kill zone mechanism and proceed with step 3.
  //    The kill zone is defined by cells sharing a house with the current strong link endpoint and EVERY previous strong
  //    link starting cell in the sequence.
  //    The X-chain candidate can be removed from the kill zone, except both endpoints ?? (other points within the
  //    x-chain can be eliminated)
  // 5. Special cases:
  //    - strong link ends in same house as root cell: sequence finishes as a "continuous loop",
  //      where all weak links in chain become effectively "strong links" and allow several eliminations
  //      (Sudoku Swami, tutorial #27, mark 35:25, mark 36:50; tutorial #42 mark 20:20)
  //    - strong link ends (and finishes) with root cell: the X-chain candidate can be established(!) in this cell.
  //      (dxSudoku #41, mark 9:45, mark 13:07, mark 13:31, mark 13:42)
  //    Long chain examples: dxSudoku #41, mark 13:20, mark 13:42.
  //
  bool found = false;

  for (int clueIndex = 0; clueIndex < SUDOKU_GRID_SIZE; clueIndex++) {
    // clear chain flags necessary?
    if (X_Chain::seedChain( RowIterator { (*grid_).data() }, clueIndex ))
      found = true;
    if (X_Chain::seedChain( ColIterator { (*grid_).data() }, clueIndex ))
      found = true;
    if (X_Chain::seedChain( BoxIterator { (*grid_).data() }, clueIndex ))
      found = true;
  }
  for (auto &cell : *grid_)
    if (cell.candidateProtected()) {
      // printf( "-------- Protected cell in *grid_ detected after x_Chain()\n" );
      assert( false );
    }

  return found;
}  // -------------------------------------------------------------------------------------------

// ==============================================================================================

void SudokuGrid::X_NodeChain::Link::Node::print( const char *prefix_format, ... ) const {

  assert( validNodeId( type_ ) );

  va_list args;
  va_start( args, prefix_format );
  vprintf( prefix_format, args );
  va_end( args );

  auto nodeSeparator = [this]() { return (isPointNode()) ? "" : (isRowGroupNode()) ? "=" : "\xba"; };
                                                           // note: '\xba' == 'two vertical lines in a single char'

  printf( "[%s%d,%d%s]", nodeSeparator(), rowNumber(cell_), colNumber(cell_), nodeSeparator());

}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::Node::printTitled( const char *prefix, const char *suffix_format, ... ) const {

  assert( validNodeId( type_ ) );

  print( "%s", prefix );

  va_list args;
  va_start( args, suffix_format );
  vprintf( suffix_format, args );
  va_end( args );

}  // -------------------------------------------------------------------------------------------

size_t SudokuGrid::X_NodeChain::Link::length( void ) const {
  // return links in chain quantity
  size_t count = 0;

  for (const Link *link = this; link != nullptr; link = link->prev_)
    count += 1;

  return count;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::printChainEngine( X_NodeChain *chain, const char *format, va_list &args ) const {

  for (const Link *link = prev_; link != nullptr; link = link->prev_->prev_) {
    link->node.print( " <-- " );
    link->prev_->node.print( " <== " );
  }
  vprintf (format, args);

}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::printChain( X_NodeChain *chain, const char *format, ... ) const {
  // print currently chained links finished by node in reversed order.
  // current link object is expected to contain recently chained weakLink ending node.

  printf( "Gr-nodes X-chain:%zd, clue:%d: ", length(), chain->clueIndex_ + 1 );
  this->node.print( "" );              // show last _chained_ node

  va_list args;
  va_start( args, format );
  printChainEngine( chain, format, args );
  va_end( args );

}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::printChain( X_NodeChain *chain, const Node &node, const char *format, ... ) const {
  // node and last chained link (i.e. current object) are assumed to be a strong link.
  // print currently chained links finished by node in reversed order.
  // current link object is expected to contain recently chained weakLink ending node.

  // start with unchained node from argument list:
  node.print( "Gr-nodes X-chain:%zd, clue:%d: ", length() + 1, chain->clueIndex_ + 1);
  this->node.print( " <== " );         // show last _chained_ node

  va_list args;
  va_start( args, format );
  printChainEngine( chain, format, args );
  va_end( args );

  if (length() + 1 > 12)
    printf( "+ + + Very long X-chain found ! ! !\n" );

}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameRow( const Node *node ) const {
  assert( this != node  &&  cell_ != node->cell_ );

  assert( validNodeId( type_ )  &&  validNodeId( node->type_ ) );

  if (isColGroupNode()  ||  node->isColGroupNode())
    return false;
  return cell_->sharingSameRow( node->cell_ );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameCol( const Node *node ) const {
  assert( this != node  &&  cell_ != node->cell_ );

  assert( validNodeId( type_ )  &&  validNodeId( node->type_ ) );

  if (isRowGroupNode()  ||  node->isRowGroupNode())
    return false;
  return cell_->sharingSameCol( node->cell_ );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameBox( const Node *node ) const {
  assert( this != node  &&  cell_ != node->cell_ );

  assert( validNodeId( type_ )  &&  validNodeId( node->type_ ) );

  return cell_->sharingSameBox( node->cell_ );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameBox_sameEntityAllowed( const Node *node ) const {
  // return true, if current object and given node share same box
  // however, in contrast to Node::sharingSameBox( const Node * ) method asking for nodes with same starting cell is allowed.

  assert( this != node );
//assert( cell_ != node->cell_ );  --> such situation is intentionally allowed (eg. for perpendicular nodes)

  assert( validNodeId( type_ )  &&  validNodeId( node->type_ ) );

  return cell_->sharingSameBox_sameEntityAllowed( node->cell_ );
}

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameHouse( const Node *node ) const {
  assert( this != node  &&  cell_ != node->cell_ );

  assert( validNodeId( type_ )  &&  validNodeId( node->type_ ) );

  // this and node are different objects and belong to the same house (row/col/box)
  return (sharingSameRow( node )  ||  sharingSameCol( node )  ||  sharingSameBox( node ));
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::sharingSameHouse( const Cell *cell ) const {

  assert( isPartOfNode( cell ) == false );

  switch (type_) {
  case NodeId::point:
    return (cell->rowOfCellsStart_ == cell_->rowOfCellsStart_  ||
            cell->colOfCellsStart_ == cell_->colOfCellsStart_  ||
            cell->boxOfCellsStart_ == cell_->boxOfCellsStart_);
  case NodeId::rowGroup:
    return (cell->rowOfCellsStart_ == cell_->rowOfCellsStart_  ||  cell->boxOfCellsStart_ == cell_->boxOfCellsStart_);
  case NodeId::colGroup:
    return (cell->colOfCellsStart_ == cell_->colOfCellsStart_  ||  cell->boxOfCellsStart_ == cell_->boxOfCellsStart_);
  default:
    assert( 0 );
    return false;
  }
}  // -------------------------------------------------------------------------------------------

const SudokuGrid::Cell *SudokuGrid::X_NodeChain::Link::pointFromGroup( const Cell *groupStart, bool rowType,
                                                                       int clueIndex, const Cell *cellToIgnore ) {
  auto &groupIterator = (rowType) ? (HouseIterator &) RowIterator { const_cast<Cell *>( groupStart ) }
                                  : (HouseIterator &) ColIterator { const_cast<Cell *>( groupStart ) };
#if 1
  for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
    const Cell *cell = &groupIterator[i];
    if (cellToIgnore != nullptr  &&  cell == cellToIgnore)
      continue;
    if (cell->candidatePossible( clueIndex ))
      return cell;
  }
#else
  if (cellToIgnore != nullptr)
    for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
      const Cell *cell = &groupIterator[i];
      if (cell != cellToIgnore  &&  cell->candidatePossible( clueIndex ))
        return cell;
    }
  else
    for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
      const Cell *cell = &groupIterator[i];
      if (cell->candidatePossible( clueIndex ))
        return cell;
    }
#endif
  assert( 0 );               // should never happen because given cells at groupStart have to contain candidate cell(s)
  return nullptr;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::diveDeeper( X_NodeChain *chain, const Node &foundNode ) const {
  // try to continue the chain with a weak link, return elimination status

  Link link { foundNode.cell_, foundNode.type_, this };    // append found node to the chain
  link.node.candidateProtection( true );

  bool found = link.addWeakLink( chain );
  link.node.candidateProtection( false );

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::chainClosedKillZone( X_NodeChain *chain, const Node &loopStartChainedNode ) const {
  // return true, if an elimination in any killZone was successful
  // killZones: every chained weakLink (traversed until 'loopStartChainedNode') defines a killZone
  // current link object must be last chained weakLink ending node.
  // loopStartChainedNode must be a chained weakLink starting node.
  bool found = false;

  //loopStartChainedNode.print( "loopStartChainedNode: " );
  //printf( " -- used as chained links backwards scan limit in chainClosedKillZone()\n" );

  for (const Link *weakLinkEnd = this; ; ) {
    auto weakLinkStart = weakLinkEnd->prev_;
    assert( weakLinkStart != nullptr  &&  weakLinkStart->prev_ != nullptr );
    //weakLinkEnd->node.print( "chainClosedKillZone():killZone " );
    //weakLinkStart->node.print( " <-- " ), printf( "\n" );

    if (killZone( chain, weakLinkEnd->node, weakLinkStart->node ))
      found = true;
    if (&weakLinkStart->node == &loopStartChainedNode)
      break;
    weakLinkEnd = weakLinkStart->prev_;
  }

  return found;
}  // -------------------------------------------------------------------------------------------

#if 0
bool SudokuGrid::X_NodeChain::Link::chainClosedKillZone( X_NodeChain *chain,
                                                         const Node &loopStartChainedNode, const Node &endNode ) const {
  // return true, if an elimination in any killZone was successful
  // apply killZones: [endNode, loopStartChainedNode] + all chained weakLinks ( traversed until loopStartChainedNode )
  bool found = false;

  printChain( chain, endNode, " chainClosedKillZone()\n" );

  loopStartChainedNode.print( "at first: loop closing weakLink: " );
  endNode.print( " <-- " ), printf( "\n" );
  if (killZone( chain, endNode, loopStartChainedNode ))    // weakLink which effectively closes the chain
    found = true;
  if (chainClosedKillZone( chain, loopStartChainedNode ))  // chained weakLinks (until loopStartChainedNode) are killZones too
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::X_NodeChain::Link::killZone( X_NodeChain *chain, const Node &nodeA, const Node &nodeB,
                                              HouseIterator &&houseIterator ) {
  // return true, if any elimination performed in along given 'houseIterator' was successful.
  // given 'endNode' and 'startNode' cells itself are excluded from elimination process.
  assert( (houseIterator.isRowType() && nodeA.isColGroupNode() == false)  ||
          (houseIterator.isColType() && nodeA.isRowGroupNode() == false)  ||  houseIterator.isBoxType() );
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {
    Cell *cell = &houseIterator[i];
    if (cell->candidateImpossible( clueIndex )  ||  nodeA.isPartOfNode( cell )  ||  nodeB.isPartOfNode( cell ))
      continue;
    if (nodeB.sharingSameHouse( cell )) {
      assert( cell->unsolved() );      // dziwny assert (po co by³ w wersji x-chain operuj¹cej bez 'nodes' ?)
      // intentionally ignore cell->candidateProtected() status, because it's used by x_chain() for chain trackking purpose
      cell->candidateCrossOut( clueIndex );

      // todo! w przysz³oœci: zanotowaæ do wykreœlenia w strukturze *chain, zamiast natychmiastowego wykreœlania

      nodeA.print( "chainKillZone( clue:%d zone:[", clueIndex + 1 );
    #if 1
      nodeB.printTitled( "::", "]) crossed out @[row:%d, col:%d]\n", rowNumber( cell ), colNumber( cell ) );
    #else
      nodeB.print( "::" );
      printf( "]) crossed out @[row:%d, col:%d]\n", rowNumber( cell ), colNumber( cell ) );
    #endif
      found = true;
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::killZone( X_NodeChain *chain, const Node &nodeA, const Node &nodeB ) {
  // look for intersections between both given nodes in all three possible houses
  bool found = false;
  
  assert( validNodeId( nodeA.type_ ) );

  if (nodeA.isColGroupNode() == false)
    if (killZone( chain, nodeA, nodeB, RowIterator { nodeA.cell_->rowOfCellsStart_ } ))
      found = true;
  if (nodeA.isRowGroupNode() == false)
    if (killZone( chain, nodeA, nodeB, ColIterator { nodeA.cell_->colOfCellsStart_ } ))
      found = true;
  if (true)
    if (killZone( chain, nodeA, nodeB, BoxIterator { nodeA.cell_->boxOfCellsStart_ } ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::chainKillZone( X_NodeChain *chain, const Node &endNode ) const {
  // apply kill zone between 'endNode' and EVERY strong link starting node chained until now.
  // current link object must be a chained weak link ending node and builds with 'endNode' current strong link.
  bool found = false;
  auto link = this;

  do {
    assert( link->prev_ != nullptr );
    link = link->prev_->prev_;                   // previous strong link start node
    if (killZone( chain, endNode, link->node ))  // kill zone is defined by endNode and each chained strong link starting node
      found = true;
  } while (link->prev_ != nullptr);

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::Node::isPartOfNode( const Cell *cell ) const {
  // return true, if given cell belongs to current node

  assert( validNodeId( type_ ) );

  switch (type_) {
  case NodeId::point:
    return cell_ == cell;
  case NodeId::rowGroup:
    return (cell_->boxOfCellsStart_ == cell->boxOfCellsStart_  &&  cell_->rowOfCellsStart_ == cell->rowOfCellsStart_);
  case NodeId::colGroup:
    return (cell_->boxOfCellsStart_ == cell->boxOfCellsStart_  &&  cell_->colOfCellsStart_ == cell->colOfCellsStart_);
  default:
    assert( 0 );
    return false;
  }
}  // -------------------------------------------------------------------------------------------

const SudokuGrid::X_NodeChain::Link::Node * SudokuGrid::X_NodeChain::Link::getChainedNode( const Cell *cell ) const {
  // return chained node containing given cell
  assert( cell->candidateProtected() );

  for (auto link = this; link != nullptr; link = link->prev_)
    if (link->node.isPartOfNode( cell ))
      return &link->node;

  assert( 0 );  // cell seems to be NOT linked
  return nullptr;
}  // -------------------------------------------------------------------------------------------

const SudokuGrid::X_NodeChain::Link *
SudokuGrid::X_NodeChain::Link::chainEffectivelyClosed( const Node &newStrongLinkEnd ) const {
  // find deepest chained strongLink starting node which shares a house with 'newStrongLinkEnd'
  // and return chained ending node of this found strongLink (which is also next weakLink starting node)
  // 
  // current link object must be a chained weakLink ending node and builds with 'newStrongLinkEnd' current strong link.
  // current chain must contain at least one weakLink/strongLink set.

  assert( prev_ != nullptr  &&  prev_->prev_ != nullptr );           // at least one chained weakLink/strongLink set

  const Link *weakLinkStart = nullptr;
  const Link *closingStrongLink = nullptr;

  for (const Link *strongLinkStart = this; strongLinkStart->prev_ != nullptr; ) {
    weakLinkStart   = strongLinkStart->prev_;
    assert( weakLinkStart->prev_ != nullptr );
    strongLinkStart = weakLinkStart->prev_;                          // 2 links back to previous strongLink start
    if (newStrongLinkEnd.sharingSameHouse( &strongLinkStart->node ))
      closingStrongLink = weakLinkStart;
  }
  return closingStrongLink;  // link [end, start] == [closingStrongLink, closingStrongLink->prev_]
#if 0
  const Link *chainedNodeInSameHouse = nullptr;

  for (auto link = this; link->prev_ != nullptr; ) {
    assert( link->prev_->prev_ != nullptr );
    link = link->prev_->prev_;
    if (strongLinkEnd.sharingSameHouse( &link->node ))               // look at each chained strongLink starting node
      chainedNodeInSameHouse = link;
  }
  return chainedNodeInSameHouse;
#endif
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::isStrongLinkStartInChain( const Node &givenNode ) const {
  // return true, if 'givenNode' is a chained strongLink node.
  // current link object must be a chained weakLink ending node and builds with 'givenNode' current strongLink.
  // 'givenNode' must belong to chain starting at current link object.

  assert( givenNode.isPointNode()  &&  givenNode.cell_->candidateProtected() );

  for (const Link *link = prev_; link != nullptr; ) {
    if (&link->node == &givenNode)
      return false;                              // found node is a strong link ENDING side

    assert( link->prev_ != nullptr );
    link = link->prev_;

    if (&link->node == &givenNode)               // found node is a strong link STARTING side
      return true;

    link = link->prev_;
  }
  assert( false );                               // 'node' not in chain??!
  return false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::usualEliminationAndDiveDeeper( X_NodeChain *chain, Node &foundNode ) const {
  // perform usual elimination according to found fresh strongLink ('last chained node' ==> 'foundNode').
  // killZones: 'foundNode' + each chained strongLink starting node.
  // if foundNode shares a house with any of chained strongLink starting nodes,
  // then the chain has a closed loop and each weakLink in this loop defines additional killZone.
  // current link object must be a chained weakLink ending node and builds with 'foundNode' current strongLink.
  bool found = false;
#if 0
  if (length() + 1 > 12)
    printChain( chain, foundNode,
                "\nusual elimination with fresh strongLink (%sNode ending)\n", (foundNode.isPointNode()) ? "Point" : "Group" );
#endif
  if (chainKillZone( chain, foundNode )) {       // usual killZones: 'foundNode' + each chained strongLink starting node
    printChain( chain, foundNode, "\nusualElimination(), last strongLink ==> %sNode (case #1)\n",
                (foundNode.isPointNode()) ? "point" : "group" );
    found = true;                                // <-- most commonly observed elimination path in x-chain!
  }
  auto closingStrongLink = chainEffectivelyClosed( foundNode );   // [end, start] == [closingStrongLink, closingStrongLink->prev_]
  if (closingStrongLink != nullptr) {
  //foundNode.print( "additionally: foundNode " );
  //closingStrongLink->prev_->node.print( " shares a house with chained strongLink starting node " );
  //printf( " =>> CLOSED LOOP!\n" );
    // killZones: 'foundNode'+'closingStrongLink->prev_->node'
    //            _AND_ all chained weakLinks traversed backwards until closingStrongLink->node
    if (killZone( chain, foundNode, closingStrongLink->prev_->node )) {
      printChain( chain, foundNode,
                  "\nusualElimination(), chainEffectivelyClosed, last strongLink ==> %sNode (case #2)\n",
                  (foundNode.isPointNode()) ? "point" : "group" );
      /* should never happen, because the killZone(chain, foundNode, loopStartLink->node) */
      /* has been already done by above chainKillZone( chain, foundNode ); */
      assert( 0 );
      found = true;
    }
    if (chainClosedKillZone( chain, closingStrongLink->node )) {  // all chained weakLinks traversed backwards until...
      printChain( chain, foundNode,
                  "\nusualElimination(), chainEffectivelyClosed, last strongLink ==> %sNode (case #3)\n",
                  (foundNode.isPointNode()) ? "point" : "group" );
      found = true;                              // <-- rare, but indeed observed elimination!
    }
  }
  if (diveDeeper( chain, foundNode ))
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::specialElimination( X_NodeChain *chain, Node &foundNode ) const {
  // special elimination occurs when 'foundNode' belongs already to chained nodes.
  // If 'foundNode' is a chained strongLink starting node, then:
  //   if 'foundNode.isPointNode()', then 'foundNode.cell_' can be immediately solved,
  //   if 'foundNode.isGroupNode()', then all candidates in house of 'foundNode' (except 'foundNode' cells) can be eliminated.
  // If 'foundNode' is a chained strongLink ending node, then each already chained weakLink is a killZone!
  // current link object must be a chained weakLink ending node and builds with 'foundNode' current strongLink.
  // further chain building (== 'diveDepper()') is NOT possible, because 'foundNode' cell(s) is/are already chained.
  bool found = false;

  if (foundNode.isGroupNode()) {
    printChain( chain, foundNode, "\n + + + Found an already completely protected group node. Is it a closed loop?\n" );
    // is such scenario possible at all???
    printf( "+ + + What can be done? More investigation required!\n" );

    assert( 0 );  /* to be implemented, (similarly to foundNode.isPointNode(), but steps need to be adapted) */
    // jeœli wsród ogniw w ³añcuchu jest ju¿ foundNode (groupNode w tym samym kierunku i tymi samymi komórkami)
    //   jeœli foundNode jest chained strongLink starting node,
    //     to nie mo¿na wskazaæ komórki w 'foundNode' do wpisania gotowej wartoœci, ale mo¿na poza 'foundNode' wyeliminowaæ
    //     pozosta³ych kandydatów w tym samym house (czyli w box i w wierszu (jeœli foundNode.isRowGroupNode()) b¹dŸ
    //     w kolumnie (jeœli foundNode.isCOlGroupNode())
    //   jeœli foundNode jest chained strongLink ending node
    //     to mamy closed chain, czyli ka¿dy weakLink w ³añcuchu definiuje kilZone!
  } else {
    const Node *chainedNode = getChainedNode( foundNode.cell_ );
    assert( chainedNode != nullptr );
    if (chainedNode->isPointNode()) {
      assert( foundNode.cell_ == chainedNode->cell_ );
      if (isStrongLinkStartInChain( *chainedNode )) {
        // found node (strongLink ending node) is also a chained strongLink starting node == nice loop case!
        printChain( chain, foundNode, "\nfound POINT node is already a chained strongLink start node ==> cell SOLVED!\n" );
        establishClueAndCrossItOutInRowColBox( const_cast<Cell *>(chainedNode->cell_), chain->getClueIndex() );
        found = true;
      } else {                                             // closed loop detected!
        //printChain( chain, foundNode, "\nfound POINT node already in chain --> closed loop!\n" );
        if (chainKillZone( chain, *chainedNode )) {        // killZones: 'chainedNode' + each chained strongLink starting node
          printChain( chain, foundNode, "\n" );
          foundNode.printTitled( "found pointNode: ", "already in chain --> (partially) closed loop!(CASE #1)\n" );
          found = true;      // indeed observed elimination!
        }
        if (chainClosedKillZone( chain, *chainedNode )) {  // killZones: chained weakLinks traversed backwards until chainedNode
          printChain( chain, foundNode, "\n" );
          foundNode.printTitled( "found pointNode: ", "already in chain --> (partially) closed loop!(CASE #2)\n" );
          found = true;      // indeed observed elimination!
        }
      }
    }
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::eliminate( X_NodeChain *chain, Node &foundNode ) const {
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  assert( validNodeId( node.type_ )  &&  validNodeId( foundNode.type_ ) );

  if (foundNode.isGroupNode()) {
    auto &groupIt = (foundNode.isRowGroupNode()) ? (HouseIterator &) RowIterator { const_cast<Cell *>( foundNode.cell_ ) }
                                                 : (HouseIterator &) ColIterator { const_cast<Cell *>( foundNode.cell_ ) };
    if (candidateGroupProtected( groupIt )) {              // group node already in chain?
      // todo!:  assert( foundNode faktycznie znajduje siê ju¿ w ³añcuchu );
      if (specialElimination( chain, foundNode ))          // closed loop, perhaps nice loop case!
        found = true;
    } else if (candidateGroupUnprotected( groupIt, clueIndex )) {
      if (usualEliminationAndDiveDeeper( chain, foundNode ))
        found = true;
    } else {
      // dead end: found group node cells are already both: protected _AND_ unprotected
      //printChain( chain, foundNode, "\nfound potential group node cells are already PARTIALLY chained (== protected)\n" );
      //printf( "" ); // investigation required == BREAK POINT!
    }
  } else if (foundNode.cell_->candidateProtected()) {      // foundNode (of point type) already in x-chain?
    if (specialElimination( chain, foundNode ))
      found = true;          // indeed observed elimination!
  } else {
    if (usualEliminationAndDiveDeeper( chain, foundNode ))
      found = true;
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addStrongLinkWithinLine( X_NodeChain *chain, HouseIterator &lineIterator ) const {
  const int clueIndex = chain->getClueIndex();
  assert( lineIterator.isLineType() );
  assert( node.isPointNode() ||
          node.isRowGroupNode() && lineIterator.isRowType() ||
          node.isColGroupNode() && lineIterator.isColType() );

  const Cell *group = node.isPointNode() ? groupNodeStartCell( node.cell_, lineIterator.isRowType() ) : node.cell_;
  struct FreshGroup {
    Node node;
    int  candidateCount;
  } freshGroup { { NodeId::point, nullptr }, 0 };

  int houseCandidateCount = 0;
  for (size_t i = 0; i < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; i++, lineIterator.nextGroup()) {
    int traversedGroupCandidateCount = candidateCountInGroup( lineIterator, clueIndex );
    if (traversedGroupCandidateCount == 0)                 // dull, empty group node
      continue;
    houseCandidateCount += traversedGroupCandidateCount;   // count all candidates in house
    if (group == &lineIterator[0]) {
      if (node.isPointNode() && traversedGroupCandidateCount == 2)
        freshGroup = { { NodeId::point, pointFromGroup( group, lineIterator.isRowType(), clueIndex, node.cell_ ) }, 1 };
    } else if (traversedGroupCandidateCount == 1)
      freshGroup = { { NodeId::point, pointFromGroup( &lineIterator[0], lineIterator.isRowType(), clueIndex ) }, 1 };
    else {
      static_assert( NodeId::rowGroup == static_cast<NodeId>( HouseId::row )  &&         // required by following
                     NodeId::colGroup == static_cast<NodeId>( HouseId::col ) );          //  node type argument short cut:
      freshGroup = { { static_cast<NodeId>( lineIterator.type() ), &lineIterator[0] }, traversedGroupCandidateCount };
    }
  }
  if (freshGroup.candidateCount == 0)
    return false;                                          // strong link not possible, no fresh group found
  if (houseCandidateCount != candidateCountInNode( node, clueIndex ) + freshGroup.candidateCount)
    return false;                                          // strong link not possible, more than 2 groups contain candidates

  return eliminate( chain, freshGroup.node );              // try to eliminate, append fresh STRONG LINK, dive deeper
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addStrongLinkWithinBox( X_NodeChain *chain, HouseIterator &&boxGroupIterator ) const {
  const int clueIndex = chain->getClueIndex();
  assert( boxGroupIterator.isLineType() );
  assert( node.isPointNode() || node.isRowGroupNode() || node.isColGroupNode() );        // i.e. node may be of ANY type

  struct FreshGroup {
    Node node;
    int  candidateCount;
  } freshGroup { { NodeId::point, nullptr }, 0 };

  auto potentialNode = [clueIndex]( HouseIterator &groupIterator,
                                    int candidateCount, const Cell *toIgnore = nullptr ) -> struct FreshGroup {
    NodeId     nodeType;
    const Cell *nodeCell;

    if (candidateCount > 1) {
      assert( groupIterator.isLineType() );
      static_assert(NodeId::rowGroup == static_cast<NodeId>(HouseId::row)  &&            // required by following
                    NodeId::colGroup == static_cast<NodeId>(HouseId::col));              //  node type argument short cut:
      nodeType = static_cast<NodeId>( groupIterator.type() );
      nodeCell = &groupIterator[0];
    } else {
      nodeType = NodeId::point;
      nodeCell = pointFromGroup( &groupIterator[0], groupIterator.isRowType(), clueIndex, toIgnore );
    }
    return { { nodeType, nodeCell }, candidateCount };
  };
  const Cell *group = node.isPointNode() ? groupNodeStartCell( node.cell_, boxGroupIterator.isRowType() ) : node.cell_;

  int houseCandidateCount = 0;
  for (size_t i = 0; i < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; i++, boxGroupIterator.nextBoxGroup()) {
    int traversedGroupCandidateCount = candidateCountInGroup( boxGroupIterator, clueIndex );
    if (traversedGroupCandidateCount == 0)                 // dull, empty group node
      continue;
    houseCandidateCount += traversedGroupCandidateCount;   // count all candidates in box
    if (node.isGroupNode()) {
      if (node.isParallel( boxGroupIterator.type() )) {
        if (group != &boxGroupIterator[0])
          freshGroup = potentialNode( boxGroupIterator, traversedGroupCandidateCount );  // catch potentially interesting group
      } else if (node.intersectionOfNodes( boxGroupIterator )->candidateImpossible( clueIndex ))
        freshGroup = potentialNode( boxGroupIterator, traversedGroupCandidateCount );    // catch potentially interesting group
      else if (traversedGroupCandidateCount == 2)
        freshGroup = potentialNode( boxGroupIterator, 1, node.intersectionOfNodes( boxGroupIterator) );
    } else if (group == &boxGroupIterator[0]) {            // node.cell_ belongs to currently traversed group
      if (traversedGroupCandidateCount == 2)               // two lonely candidates in box?
        freshGroup = potentialNode( boxGroupIterator, 1, node.cell_ );
    } else                                                 // node.cell_ DOES NOT belong to currently traversed group
      freshGroup = potentialNode( boxGroupIterator, traversedGroupCandidateCount );      // catch potentially interesting group
  }

  if (freshGroup.candidateCount == 0)
    return false;                                          // strong link not possible, no fresh group found
  if (houseCandidateCount != candidateCountInNode( node, clueIndex ) + freshGroup.candidateCount)
    return false;                                          // strong link not possible, more than 2 groups contain candidates

  return eliminate( chain, freshGroup.node );              // try to eliminate, append fresh STRONG LINK, dive deeper
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addStrongLink( X_NodeChain *chain, HouseIterator &houseIterator ) const {
  bool found = false;

  assert( validHouse( houseIterator.type() ) );
  assert( node.isPointNode() ||
          node.isRowGroupNode() && houseIterator.isColType() == false ||
          node.isColGroupNode() && houseIterator.isRowType() == false );

  if (houseIterator.isLineType()) {
    assert( node.isPointNode() ||
            node.isRowGroupNode() && houseIterator.isRowType() ||
            node.isColGroupNode() && houseIterator.isColType() );
    if (addStrongLinkWithinLine( chain, houseIterator ))                                  // search in given row/column
      found = true;
  } else if (node.isPointNode() && candidateCountInHouse( houseIterator, chain->getClueIndex() ) == 2) {
    assert( houseIterator.isBoxType() );
    assert( node.cell_->occupyingDifferentBox( prev_->node.cell_ ) );
    if (addStrongLinkWithinBox( chain, RowIterator { &houseIterator[0] } ))               // search in all row groups only
      found = true;
    // strong link searching in column groups would find same result like above row group search
  } else {
    assert( houseIterator.isBoxType() );
    assert( node.isPointNode()    && node.cell_->occupyingDifferentBox( prev_->node.cell_ ) ||
            node.isRowGroupNode() && node.cell_->sharingSameRow( prev_->node.cell_ ) ||
            node.isColGroupNode() && node.cell_->sharingSameCol( prev_->node.cell_ ) );
    if (addStrongLinkWithinBox( chain, RowIterator { &houseIterator[0] } ))               // search in all row groups
      found = true;
    if (addStrongLinkWithinBox( chain, ColIterator { &houseIterator[0] } ))               // search in all col groups
      found = true;
  }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addStrongLink( X_NodeChain *chain ) const {
  // append any possible "strong link", apply kill zone and go recursively further 
  bool found = false;

  // try to find each possible strong link outgoing from "this->node"
  // if cross out successful in kill zone, found = true;
  // add found "strong link" node (cell or group) and dive deeper

  assert( validNodeId( node.type_ ) );

  // make sure that nodes of last chained weakLink are not same:
  assert( node.isPointNode() && prev_->node.isPointNode() && node.cell_ != prev_->node.cell_  ||
          node.isPointNode() && prev_->node.isGroupNode() && prev_->node.isPartOfNode( node.cell_ ) == false  ||
          node.isGroupNode() && prev_->node.isPointNode() && node.isPartOfNode( prev_->node.cell_ ) == false  ||
          node.isGroupNode() && prev_->node.isGroupNode() &&
          (node.isParallel( prev_->node ) && node.cell_ != prev_->node.cell_ ||
           node.isPerpendicular( prev_->node ) && (node.occupyingDifferentBox_sameEntityAllowed( &prev_->node ) ||
                                                   node.intersectionOfNodes( prev_->node )->candidateImpossible( chain->getClueIndex() ))) );
  if (node.isPointNode()) {
    if (node.occupyingDifferentRow( &prev_->node ))
      if (addStrongLink( chain, RowIterator { node.cell_->rowOfCellsStart_ } ))
        found = true;
    if (node.occupyingDifferentCol( &prev_->node ))
      if (addStrongLink( chain, ColIterator { node.cell_->colOfCellsStart_ } ))
        found = true;
    if (node.occupyingDifferentBox( &prev_->node ))
      if (addStrongLink( chain, BoxIterator { node.cell_->boxOfCellsStart_ } ))
        found = true;
  } else if (node.sharingSameBox_sameEntityAllowed( &prev_->node )) {
    if (addStrongLink( chain, (node.isRowGroupNode()) ? (HouseIterator &) RowIterator { node.cell_->rowOfCellsStart_ }
                                                      : (HouseIterator &) ColIterator { node.cell_->colOfCellsStart_ }))
      found = true;
  } else if (node.isRowGroupNode() && node.cell_->sharingSameRow( prev_->node.cell_ )  ||
             node.isColGroupNode() && node.cell_->sharingSameCol( prev_->node.cell_ )) {
    if (addStrongLink( chain, BoxIterator { node.cell_->boxOfCellsStart_ } ))
      found = true;
  } else
    assert( 0 );                                 // last link seems to be invalid

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addGroupedNodeWithinLineWeakLink( X_NodeChain *chain, HouseIterator &lineIterator ) const {
  bool found = false;
  const int clueIndex = chain->getClueIndex();
  assert( lineIterator.isLineType() );

  for (size_t box = 0; box < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; box++, lineIterator.nextGroup())
    if (candidateGroupUnprotected( lineIterator, clueIndex ) && candidateCountInGroup( lineIterator, clueIndex ) > 1) {
      // append found fresh grouped node (also not being current node itself), then go recursively deeper
      static_assert(NodeId::rowGroup == static_cast<NodeId>(HouseId::row)  &&  // required by following Link
                    NodeId::colGroup == static_cast<NodeId>(HouseId::col));    //  constructor's argument short cut:
      Link link { /* groupStart: */ &lineIterator[0], static_cast<NodeId>(lineIterator.type()), this};
      link.node.candidateProtection( true );
    //link.node.print( "weakLink --> " ), printf( "\n" );
      if (link.addStrongLink( chain ))
        found = true;
      link.node.candidateProtection( false );
    }
  //lineIterator.reset();

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addGroupedNodeWithinBoxWeakLink( X_NodeChain *chain, HouseIterator &&boxGroupIterator ) const {
  // return true, if underlying code performs successful elimination
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  assert( boxGroupIterator.isLineType() );
  bool (*fineAsWeakLink)(const Node & node, HouseIterator & boxGroupIterator, int clueIndex);
  auto freshNodeUprotectedAndInDifferentRow = [](const Node &node, HouseIterator &boxGroupIterator, int clueIndex) {
    return candidateGroupUnprotected( boxGroupIterator, clueIndex )  &&
           node.cell_->rowOfCellsStart_ != boxGroupIterator[0].rowOfCellsStart_;
  };
  auto freshNodeUprotectedAndInDifferentCol = []( const Node &node, HouseIterator &boxGroupIterator, int clueIndex ) {
    return candidateGroupUnprotected( boxGroupIterator, clueIndex )  &&
           node.cell_->colOfCellsStart_ != boxGroupIterator[0].colOfCellsStart_;
  };
  auto freshNodeUnprotectedAndCollisionFree = []( const Node &node, HouseIterator &boxGroupIterator, int clueIndex ) {
    const Cell *nodeIntersection = node.isRowGroupNode() ? intersectionOfCells( node.cell_, &boxGroupIterator[0] )
                                                         : intersectionOfCells( &boxGroupIterator[0], node.cell_ );
    if (nodeIntersection->candidatePossible( clueIndex ))
      return false;                              // current node and fresh node intersection cannot contain a candidate
    for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
      const Cell *freshNodeCell = &boxGroupIterator[i];
      if (freshNodeCell != nodeIntersection && freshNodeCell->candidateProtected())
        return false;                            // cell seems to be already in chain
    }
    return true;         // fresh node cells are unprotected AND node intersection does NOT contain current candidate
  };

  if (boxGroupIterator.isRowType() && (node.isPointNode() || node.isRowGroupNode()))
    fineAsWeakLink = freshNodeUprotectedAndInDifferentRow;
  else if (boxGroupIterator.isColType() && (node.isPointNode() || node.isColGroupNode()))
    fineAsWeakLink = freshNodeUprotectedAndInDifferentCol;
  else {                                         // fresh node perpendicularly to current node
    assert( node.isGroupNode() && node.isPerpendicular( boxGroupIterator.type() ) );
    fineAsWeakLink = freshNodeUnprotectedAndCollisionFree;
  }

  for (size_t groupNo = 0; groupNo < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; groupNo++, boxGroupIterator.nextBoxGroup()) {
    const Cell *groupStart = &boxGroupIterator[0];
    if (fineAsWeakLink( node, boxGroupIterator, clueIndex )  &&
        candidateCountInGroup( boxGroupIterator, clueIndex ) > 1) {
    // append found fresh grouped node (also not being current node itself), then go recursively deeper
      static_assert(NodeId::rowGroup == static_cast<NodeId>(HouseId::row)  &&          // required by following Link
                    NodeId::colGroup == static_cast<NodeId>(HouseId::col));            //  constructor's argument short cut:
      Link link { groupStart, static_cast<NodeId>(boxGroupIterator.type()), this };
      link.node.candidateProtection( true );
    //link.node.print( "weakLink --> " ), printf("\n");
      if (link.addStrongLink( chain ))
        found = true;
      link.node.candidateProtection( false );
    }
  }

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addWeakLink( X_NodeChain *chain, HouseIterator &houseIterator ) const {
  // append any possible "weak link", try with all potential cells in range of given houseIterator
  bool found = false;
  const int clueIndex = chain->getClueIndex();

  assert( validHouse( houseIterator.type() ) );
  assert( node.isPointNode()  ||
          node.isRowGroupNode() && houseIterator.isColType() == false  ||
          node.isColGroupNode() && houseIterator.isRowType() == false );

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++) {          // weak link made to a single cell (node.type_ == NodeId::point)
    Cell *cell = &houseIterator[i];
    if (cell->candidatePossible( clueIndex ) && cell->candidateUnprotected()) {
      Link link { cell, NodeId::point, this };             // append new link to the chain and dive deeper
      cell->candidateProtection( true );
      //printf( "weakLink --> [%d,%d]\n", rowNumber( cell ), colNumber( cell ) );
      if (link.addStrongLink( chain ))
        found = true;
      cell->candidateProtection( false );
    }
  }
  if (houseIterator.isLineType()) {                        // weak link made to a group node:
    assert( node.isPointNode()  ||
            node.isGroupNode() && node.cell_->sharingSameBox( prev_->node.cell_ ) );
    if (addGroupedNodeWithinLineWeakLink( chain, houseIterator ))
      found = true;
  } else {                                                 // weak link (from any node) to a group node within same box
    assert( houseIterator.isBoxType() );
    assert( node.isPointNode()    && node.cell_->occupyingDifferentBox( prev_->node.cell_ ) ||
            node.isRowGroupNode() && node.cell_->sharingSameRow( prev_->node.cell_ ) ||
            node.isColGroupNode() && node.cell_->sharingSameCol( prev_->node.cell_ ) );
    // reject weak links in a box between two pependicular group nodes, if these have common candidate on own junction,
    // accept remaining possibilities (parallel and perpendicular group nodes within current box)
    if (addGroupedNodeWithinBoxWeakLink( chain, RowIterator { &houseIterator[0] } ))     // search across all row groups
      found = true;
    if (addGroupedNodeWithinBoxWeakLink( chain, ColIterator { &houseIterator[0] } ))     // search across all col groups
      found = true;
  }
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::Link::addWeakLink( X_NodeChain *chain ) const {
  // append any possible "weak link", skip searching in the house occupied by last strong link
  bool found = false;

  //printf( "addWeakLink( clue:%d ) Either-Or link: @[%d,%d ==> %d,%d]\n", chain->getClueIndex() + 1,
  //        rowNumber( prev_->cell_ ), colNumber( prev_->cell_ ),
  //        rowNumber( cell_ ),        colNumber( cell_ ) );
  assert( validNodeId( node.type_ ) );

  if (node.isPointNode()) {                      // skip the house occupied by last strong link
    if (node.occupyingDifferentRow( &prev_->node ))
      if (addWeakLink( chain, RowIterator { node.cell_->rowOfCellsStart_ } ))
        found = true;
    if (node.occupyingDifferentCol( &prev_->node ))
      if (addWeakLink( chain, ColIterator { node.cell_->colOfCellsStart_ } ))
        found = true;
    if (node.occupyingDifferentBox( &prev_->node ))
      if (addWeakLink( chain, BoxIterator { node.cell_->boxOfCellsStart_ } ))
        found = true;
  } else if (node.sharingSameBox_sameEntityAllowed( &prev_->node )) {
    if (addWeakLink( chain, (node.isRowGroupNode()) ? (HouseIterator &) RowIterator { node.cell_->rowOfCellsStart_ }
                                                    : (HouseIterator &) ColIterator { node.cell_->colOfCellsStart_ } ))
      found = true;
  } else if (node.isRowGroupNode() && node.cell_->sharingSameRow( prev_->node.cell_ ) ||
             node.isColGroupNode() && node.cell_->sharingSameCol( prev_->node.cell_ )) {
    if (addWeakLink( chain, BoxIterator { node.cell_->boxOfCellsStart_ } ))
      found = true;
  } else
    assert( 0 );                                 // last link is invalid
  return found;
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::Node::candidateProtection( bool state ) {

  if (isGroupNode()) {
    assert( isRowGroupNode() || isColGroupNode() );
    Cell *groupStart = const_cast<Cell *>( cell_ );
    HouseIterator &groupIterator = (isRowGroupNode()) ? (HouseIterator &) RowIterator { groupStart }
                                                      : (HouseIterator &) ColIterator { groupStart };
    SudokuGrid::candidateGroupProtection( state, groupIterator );
  } else {
    assert( isPointNode() );
    cell_->candidateProtection( state );
  }
}  // -------------------------------------------------------------------------------------------

void SudokuGrid::X_NodeChain::Link::arrangeSeedLink( HouseIterator &tailGroupIt, int tailCandidateCount,
                                                     HouseIterator &headGroupIt, int headCandidateCount,
                                                     int clueIndex, const Cell *intersectionCell ) {
  // arrange seedLink to form a strongLink: [tailGroupNode <== headGroupNode]

  Link *seedLink = this;  // current object must be the first element of Link vector with 2 chained elements

  assert( seedLink[0].prev_ == &seedLink[1]  &&  seedLink[1].prev_ == nullptr );                   // linkage ok?

  seedLink[0].node = getNode( tailGroupIt, tailCandidateCount, clueIndex );
  seedLink[1].node = getNode( headGroupIt, headCandidateCount, clueIndex, intersectionCell );      // chain root

}  // -------------------------------------------------------------------------------------------

SudokuGrid::X_NodeChain::Link::Node
SudokuGrid::X_NodeChain::getNode( HouseIterator &groupIterator,
                                  int candidateCount, int clueIndex, const Cell *intersectionCell ) {
  // combine node type and pointer to: single cell / cell group (depending on type)
  // reduce groupNode to pointNode, if:
  //   ** groupNode have a "busy" intersection with another perpendicular groupNode
  //   ** groupNode candidateCount equals 2
  assert( groupIterator.isLineType() );
  assert( intersectionCell->candidatePossible( clueIndex )   && candidateCount > 1  ||
          intersectionCell->candidateImpossible( clueIndex ) && candidateCount > 0 );

  auto singleCell = [&]() -> const Cell * {
    for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
      const Cell *cell = &groupIterator[i];
      if (cell->candidatePossible( clueIndex ) && cell != intersectionCell)
        return cell;
    }
    assert( 0 );                       // should never happen
    return nullptr;
  };
  if (intersectionCell->candidatePossible( clueIndex ) && candidateCount == 2)
    return Link::Node { NodeId::point, singleCell() };     // reduce groupNode to pointNode!
  else
    return getNode( groupIterator, candidateCount, clueIndex );
}  // -------------------------------------------------------------------------------------------

SudokuGrid::X_NodeChain::Link::Node
SudokuGrid::X_NodeChain::getNode( HouseIterator &groupIterator, int candidateCount, int clueIndex ) {
  // combine node type and pointer to: single cell / cell group (depending on type)
  assert( groupIterator.isLineType() );
  assert( candidateCount > 0 );

  auto nodeType = [&]() -> NodeId {
    return (candidateCount == 1) ? NodeId::point : (groupIterator.isRowType()) ? NodeId::rowGroup : NodeId::colGroup;
  };
  auto nodeCell = [&]() -> const Cell * {
    if (candidateCount == 1)
      for (size_t i = 0; i < SUDOKU_BOX_SIZE; i++) {
        const Cell *cell = &groupIterator[i];
        if (cell->candidatePossible( clueIndex ))
          return cell;
      }
    return &groupIterator[0];
  };

  return Link::Node { nodeType(), nodeCell() };  // favor Return Value Optimization (RVO is compiler's job since C++17)
}  // -------------------------------------------------------------------------------------------

// ==== old x-chain: with fewer seedLinks than possible (in special cases) ======================

   // todo!  task #3:
   // are current addWeakLink() and addStrongLink() really able to add
   // both PERPENDICULAR and PARALLEL groupNodes within same box -- test required!

// ==============================================================================================

void SudokuGrid::X_NodeChain::printSeedLink( const Link (&seedLink)[2], const char *format, ... ) {

  va_list args;
  va_start( args, format );
  vprintf( format, args );
  va_end( args );

  seedLink[0].node.print( "" );
  seedLink[1].node.print( " <== ", "\n" );

}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::cultureChain( Link (&seedLink)[2] ) {
  // investigate in both possible directions all possible x-chains.
  // If any elimination performed, return true.
  bool found = false;

  //if (seedLink[0].node.isGroupNode() && seedLink[1].node.isGroupNode())
  //  printSeedLink( seedLink, "!!! Rarity: seedLink ( clue:%d ) with %s groupNodes: ",
  //                 getClueIndex() + 1,
  //                 (seedLink[0].node.isParallel( seedLink[1].node )) ? "__parallel__" : "PERPENDICULAR" );
  //else if (seedLink[0].node.isGroupNode() || seedLink[1].node.isGroupNode())
  //  printSeedLink( seedLink, "seedLink ( clue:%d ) with groupNode: ", getClueIndex() + 1 );

  seedLink[1].candidateProtection( true ), seedLink[0].candidateProtection( true );
  assert( seedLink[0].prev_ != nullptr  &&  seedLink[0].prev_->prev_ == nullptr  &&  seedLink[1].prev_ == nullptr );
  if (seedLink[0].addWeakLink( this ))                     // seedStrongLink head as chain start
    found = true;

  std::swap( seedLink[0].node, seedLink[1].node );         // exchange seedStrongLink nodes

  assert( seedLink[0].prev_ != nullptr  &&  seedLink[0].prev_->prev_ == nullptr  &&  seedLink[1].prev_ == nullptr );
  if (seedLink[0].addWeakLink( this ))                     // seedStrongLink tail as chain start
    found = true;

  seedLink[0].candidateProtection( false ), seedLink[1].candidateProtection( false );
  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedGroupStrongLinkPerpendicular( BoxIterator &boxIterator, int candidateCount ) {
  // look for origin strongLink (or strongLinks(!)) in given box house made of 2 perpendicular groupNodes.
  // One of both nodes may contain one candidate cell only.
  // Investigate all possible x-chains.
  // If any elimination is performed by any cultured x-chain, return true.
  // Following candidate cell constellation types lead to DOUBLED origin strongLinks:
  //  x.x
  //  x..  ==> 2 possible seedStrongLinks: rowGroupNode ==> pointNode, colGroupNode ==> pointNode
  //  ...                                  (groupNode is made of 2 candidate cells)
  const int clueIndex = getClueIndex();

  RowIterator rowGroupIt { &boxIterator[0] };              // look for perpendicular nodes
  ColIterator colGroupIt { &boxIterator[0] };

  for (size_t row = 0; row < SUDOKU_BOX_SIZE; row++, rowGroupIt.nextBoxGroup(), colGroupIt.reset())
    for (size_t col = 0; col < SUDOKU_BOX_SIZE; col++, colGroupIt.nextBoxGroup()) {
      const auto [intersectionCell,  intersectionBusy]  = intersection( &rowGroupIt[0], &colGroupIt[0], clueIndex );
      const auto [rowCandidateCount, colCandidateCount] = candidateCountInGroup( rowGroupIt, colGroupIt, clueIndex );
      //printf( "intersectionCell: @[%d:%d], node sizes: %d,%d\n",
      //        rowNumber( intersectionCell ), colNumber( intersectionCell ), rowCandidateCount, colCandidateCount );
      if (candidateCount + (intersectionBusy ? 1 : 0) > rowCandidateCount + colCandidateCount)
        continue;                                          // too many candidates in box for perpendicular nodes
      if (intersectionBusy && (rowCandidateCount == 1 || colCandidateCount == 1))
        continue;                                          // one of perpendicular nodes is empty
      if (intersectionBusy && rowCandidateCount > 2 && colCandidateCount > 2)
        continue;                                          // case useless for x-chain inference
      assert( rowCandidateCount > 0 && colCandidateCount > 0 );
      if (rowCandidateCount >= colCandidateCount)
        seedLink->arrangeSeedLink( rowGroupIt, rowCandidateCount,
                                   colGroupIt, colCandidateCount, clueIndex, intersectionCell );
      else
        seedLink->arrangeSeedLink( colGroupIt, colCandidateCount,
                                   rowGroupIt, rowCandidateCount, clueIndex, intersectionCell );
      //if (intersectionBusy == false && rowCandidateCount > 1 && colCandidateCount > 1)
      //  printf( "! ! ! seedLink: perpendicular groupNodes\n" );
      bool found = cultureChain( seedLink );     // side effect: exchanged seedLink[0].node and seedLink[1].node

      if (intersectionBusy && rowCandidateCount == 2 && colCandidateCount == 2) {        /* further seedLink ? */
        // constellation type with _TWO_ node combinations:
        //   x.x      xx.      x.x
        //   ...  or  x..  or  ..x  or  ...
        //   x..      ...      ...
        // 
        // ---> convert seedLink from [rowGroupNode <== pointNode] to [colGroupNode <== pointNode]
        assert( seedLink[0].node.isPointNode() && seedLink[1].node.isRowGroupNode() );
        seedLink->arrangeSeedLink( colGroupIt, colCandidateCount,
                                   rowGroupIt, rowCandidateCount, clueIndex, intersectionCell );
        assert( seedLink[0].node.isColGroupNode() && seedLink[1].node.isPointNode() );
        if (cultureChain( seedLink ))
          found = true;
      }
      return found;
    }
  return false;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedGroupStrongLinkParallel( HouseIterator &&groupIterator ) {
  // look for origin strongLink as parallel groupNode constellation: rowwise then columnwise
  // Each of both nodes has to contain at least 2 candidate cells, because nodes with 1 candidate cell 
  // are handled already by seedGroupStrongLinkPerpendicular() method.
  // Investigate all possible x-chains.
  // If any elimination is performed by cultured any x-chain, return true.
  assert( groupIterator.isRowType() || groupIterator.isColType() );

  size_t nodeCount = 0;
  const int clueIndex = getClueIndex();

  for (size_t i = 0; i < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; i++, groupIterator.nextBoxGroup()) {
    int groupCandidateCount = candidateCountInGroup( groupIterator, clueIndex );
    if (groupCandidateCount != 0) {
      if (nodeCount < elementsof( seedLink ))
        seedLink[nodeCount].node = getNode( groupIterator, groupCandidateCount, clueIndex );

      nodeCount += 1;
    }
  }
  if (nodeCount == elementsof( seedLink )  &&    // precisely 2 nodes traversed?
      seedLink[0].node.isGroupNode()  &&  seedLink[1].node.isGroupNode())
    return cultureChain( seedLink );
  else
    return false;                                // no valid seedStronLink found
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedGroupStrongLinkParallel( BoxIterator &boxIterator ) {
  // look for origin strongLink as parallel node constellation: rowwise then columnwise
  // Each of both nodes has to contain at least 2 candidate cells.
  // Investigate all possible x-chains.
  // If any elimination is performed by cultured any x-chain, return true.
  bool found = false;

  if (seedGroupStrongLinkParallel( RowIterator { &boxIterator[0] } ))          // parallel nodes rowwise
    found = true;
  if (seedGroupStrongLinkParallel( ColIterator { &boxIterator[0] } ))          // parallel nodes columnwise
    found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedGroupStrongLink( BoxIterator &boxIterator ) {
  // look for origin strongLink (or strongLinks(!)) in given box house, investigate all possible x-chains.
  // If any elimination is performed by cultured any x-chain, return true.
  // Following candidate cell constellation types lead to DOUBLED origin strongLinks:
  //  x.x
  //  x..  ==> 2 possible seedStrongLinks: rowGroupNode <== pointNode, colGroupNode <== pointNode
  //  ...                                  (groupNode is made of 2 candidate cells)
  // 
  //  x.x
  //  x.x  ==> 2 possible parallel seedStrongLinks: rowGroupNode <== rowGroupNode, colGroupNode <== colGroupNode
  //  ...                                           (both groupNodes are made of 2 candidate cells)
  const int candidateCount = candidateCountInHouse( boxIterator, getClueIndex() );

  if (candidateCount > SUDOKU_BOX_SIZE + SUDOKU_BOX_SIZE)            // upper limit: 2 parallel nodes
    return false;

  assert( candidateCount > 2 );                                      // at least one groupNode possible

  // perpendicular and parallel constellation exclude each other (in same box):
  return seedGroupStrongLinkPerpendicular( boxIterator, candidateCount )  ||
         seedGroupStrongLinkParallel( boxIterator );
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedGroupStrongLink( HouseIterator &lineIterator ) {
  // look for origin strongLink in given line house, investigate all possible x-chains.
  // If any elimination is performed by any cultured x-chain, return true.
  assert( lineIterator.isRowType() || lineIterator.isColType() );

  size_t nodeCount = 0;
  const int clueIndex = getClueIndex();

  auto &groupIterator = (lineIterator.isRowType()) ? (HouseIterator &) RowIterator { &lineIterator[0] }
                                                   : (HouseIterator &) ColIterator { &lineIterator[0] };
  for (size_t i = 0; i < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; i++, groupIterator.nextGroup()) {
    int groupCandidateCount = candidateCountInGroup( groupIterator, clueIndex );
    if (groupCandidateCount != 0) {
      if (nodeCount < elementsof( seedLink ))
        seedLink[nodeCount].node = getNode( groupIterator, groupCandidateCount, clueIndex );

      nodeCount += 1;
    }
  }
  if (nodeCount == elementsof( seedLink ))       // precisely 2 nodes traversed?
    return cultureChain( seedLink );
  else
    return false;                                // no valid seedStronLink found
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedStrongLink( HouseIterator &houseIterator ) {
  size_t candidateCountInHouse = 0;
  const int clueIndex = getClueIndex();

  for (size_t i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (houseIterator[i].candidatePossible( clueIndex )) {
      if (candidateCountInHouse < elementsof( seedLink ))
        seedLink[candidateCountInHouse].node = Link::Node { NodeId::point, &houseIterator[i] };
      candidateCountInHouse += 1;
    }

  if (candidateCountInHouse < elementsof( seedLink ))                // one candidate cell is not enough
    return false;                                                    // candidate solved or naked/hidden single

  if (candidateCountInHouse == elementsof( seedLink ))               // strongLink without any groupNodes
    return cultureChain( seedLink );

  if (houseIterator.isBoxType())
    return seedGroupStrongLink( (BoxIterator &) houseIterator );     // try a strongLink within a box with node/nodes
  else
    return seedGroupStrongLink( houseIterator );                     // try a strongLink using group node/nodes
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::X_NodeChain::seedChain( HouseIterator &&houseIterator, int clueIndex ) {
  // look for origin strongLink (or strongLinks(!)) in house, investigate all possible x-chains.
  // If any elimination performed, return true.
  bool found = false;

  X_NodeChain chain { clueIndex };

  for (size_t house = 0; house < SUDOKU_GRID_SIZE; house++, houseIterator.nextHouse())
    if (chain.seedStrongLink( houseIterator ))
      found = true;

  return found;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::x_NodeChain( void ) {
  // This is an extension of x_Chain() method, which operates additionally on group nodes instead of single cells only.
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: very advanced. Description:
  //   YouTube: "dxSudoku #41 X-Chain Puzzle Solving Technique" (= basic version without group nodes)
  //   YouTube: "dxSudoku #73 Improved X-Chain Search Algorithm" (= basic version without group nodes)
  //   YouTube: "dxSudoku #55 Types of Sudoku Links" (mark 10:30)
  //   YouTube: "dxSudoku #75 Empty Rectangle with Group Nodes and Group Links (marks 4:30, 5:35)"
  //   YouTube: "dxSudoku #77 Improved Finned X-Wing Search Algorithm" (mark 1:30, Sashimi chain may be easily extended!)
  //
  // strong link starts/ends with a group node - kill zone: all candidate cells except of node(s) starting/ending chain
  // * strong link starting with a node may end:
  //    - out of the starting box, but then in the same line as the starting node
  //    - within the SAME box (with a single cell or other node)
  //    
  // * node starting/finishing a strong link has only TWO houses: box and line of node (important for kill zone definition)
  // * all cells forming a node defining a kill zone must be protected from cross out, (cells of nodes at chain's start/end)
  // * a strong link may consist of none, one or two nodes in any order
  //
  bool found = false;

  for (int clueIndex = 0; clueIndex < SUDOKU_GRID_SIZE; clueIndex++) {
    if (X_NodeChain::seedChain( RowIterator { (*grid_).data() }, clueIndex ))
      found = true;
    if (X_NodeChain::seedChain( ColIterator { (*grid_).data() }, clueIndex ))
      found = true;
    if (X_NodeChain::seedChain( BoxIterator { (*grid_).data() }, clueIndex ))
      found = true;
  }
  for (auto &cell : *grid_)
    if (cell.candidateProtected()) {
      // printf( "-------- Protected cell in *grid_ detected after x_Chain()\n" );
      assert( false );
    }

  return found;
}  // -------------------------------------------------------------------------------------------

#if 0
bool SudokuGrid::remotePair( void ) {
  // return false, if neither clue possibility can be crossed out
  // Algorithm level: advanced. Description:
  //   YouTube: "dxSudoku #33 Remote Pair Puzzle Solving Technique"
  //   YouTube: "dxSudoku #79 Improved Remote Pair Search Algorithm"
  //   YouTube: "dxSudoku #95 Remote Pairs combined with Simple Colors"
  bool found = false;

  //
  // 1. count all cells with same pair of candidates (XY)
  // 2. minimum count of such pairs is 4
  // 3. chain these cells together: consecutive cells in the chain must be located in the same group (row/col/box)
  //    PROBLEM: there may be more than one separated chain and/or it may be branchy
  //             both cases result in complexer code to deal with such special cases
  // 4. mark which cells in chain are even and which are odd
  // 5. determine the kill zones, i.e. cells sharing the same group with any two cells from the chain
  // 6. remove XY-candidates from the kill zones, if they are in reach of an even AND an odd chain cells
  // 7. remote pair algorithm seems to be very similar to x-chain
  //

  return found;
}  // -------------------------------------------------------------------------------------------
#endif

bool SudokuGrid::validClue( Cell &cell, uchar clue ) {

  RowIterator rowOfCells { cell.rowOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (rowOfCells[i].clue_ == clue)
      return false;

  ColIterator colOfCells { cell.colOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (colOfCells[i].clue_ == clue)
      return false;

  BoxIterator boxOfCells { cell.boxOfCellsStart_ };
  for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
    if (boxOfCells[i].clue_ == clue)
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::bruteForce( size_t startIndex ) {
// return false, if failed to solve
  assert( startIndex <= elementsof( *grid_ ));

  static size_t iterationCount;
  if (startIndex == 0)
    iterationCount = 0;

  for (size_t i = startIndex; i < elementsof( *grid_ ); i++) {
    auto &cell = (*grid_)[i];
    if (cell.solved())
      continue;
    for (int k = 0; k < elementsof( cell.candidate_.clues_ ); k++)
      if (cell.candidate_.clues_[k] != 0) {
        if (validClue( cell, k + 1 )) {
          cell.setClue( k + 1 );       // store new clue
          iterationCount += 1;
          if (bruteForce( i + 1 ))     // solve remaining sudoku cells
            return true;
          else
            cell.resetClue();          // cancel, because solving was impossible
        }
      }
    return false;                      // failed, i.e. previously chosen clue(s) is/are wrong
  }

  printf( "bruteForce(): finished, after %zu recursive calls\n", iterationCount);
  return true;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::bruteForce( void ) {
  // return always false, i.e. no further steps are required. According to convention of other sudoku solving methods.

  printf( "bruteForce(): " );
  size_t initialClueCount = 0;
  for (auto &cell : *grid_)       // todo! do it with STL
    if (cell.solved())
      initialClueCount += 1;
  
  if (initialClueCount < 17) {
    printf( "Not enough initial clues (%zu), riddle cannot be solved\n", initialClueCount );
      return false;
  }

  if (bruteForce( 0 ) == false)
    printf( "bruteForce(): failed, sudoku grid seems to be invalid!\n" );

  return false;
}  // -------------------------------------------------------------------------------------------

#if 0
//
// std::hash<> function specialization required by unordered_multimap (due to customized key type)
//
template<>
struct std::hash<PossibilityFlags> {
  size_t operator()( const PossibilityFlags &key ) const {
    printf( "todo! std::hash() specialization for PossibilityFlags !!\n" );
    //return std::hash<std::string>()(key.name);
    return 0;
  }
};
#endif
#if 0
//
// std::less<> function specialization for map / multimap (due to customized key data type)
// 
// alternative: overloading operator< for customized key data type works also. However there is
//              a serious drawback: operator< would be overloaded universally for key type/class,
//              which is undesirable in client scenarios.
//
template<>
struct std::less<PossibilityFlags> {
  bool operator()( const PossibilityFlags &rhs ) const {
    printf( "todo! std::less() specialization for PossibilityFlags !!\n" );
    //return this-> ...  <  rhs. ...;
    return true;
  }
};
#endif

bool SudokuGrid::sortedBruteForce( std::multimap<int, Cell *> &grid,
                                   std::multimap<int, Cell *>::iterator it ) {
  // return false, if failed to solve

  static size_t iterationCount;
  if (it == grid.begin()) {
    for (; it != grid.end(); ++it)                         // skip solved cells
      if (it->second->unsolved())
        break;
    iterationCount = 0;
  }
  if (it == grid.end()) {
    printf( "sortedBruteForce(): finished, after %zu recursive calls\n", iterationCount);
    return true;
  }
  for (auto [k, cell] = std::pair { 0U, it++->second }; k < elementsof( cell->candidate_.clues_ ); k++)
    if (cell->candidate_.clues_[k] != 0  &&  validClue( *cell, k + 1 )) {
      cell->setClue( k + 1 );                              // store and try deeper
      iterationCount += 1;
      if (sortedBruteForce( grid, it ))                    // solve remaining sudoku cells
        return true;
      else
        cell->resetClue();                                 // cancel, solving was impossible
    }

  return false;                                            // previously chosen clue(s) must be wrong
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::sortedBruteForce( void ) {
  // return always false, i.e. no further steps are required. According to convention of other sudoku solving methods.
  // Algorithm description:
  // Sort all cells in grid so, that cells with less possibilities will be processed first.
  // This reduces the statistical probability of mistakes at begin which reduces the recursion quantity.

  std::multimap<int, Cell *> orderedGrid;

  for (auto &cell : *grid_)
  #if 1
    orderedGrid.insert( std::pair { cell.solved() ? 0 : cell.candidateCount(), &cell } );
  #else
    orderedGrid.insert( std::pair<int, Cell *>( cell.solved() ? 0 : cell.candidateCount(), &cell) );
  #endif

  if (sortedBruteForce( orderedGrid, orderedGrid.begin() ) == false)
    printf( "sortedBruteForce(): failed, sudoku grid seems to be invalid!\n" );

  return false;
}  // -------------------------------------------------------------------------------------------


void SudokuGrid::printSudokuGrid( const char *format, ... ) const {
  va_list args;
  va_start( args, format );
  vprintf( format, args );
  va_end( args );

  printf( "%c", DOUBLE_FRAME_UL_CORNER );                  // upper frame
  for (int i = 1; i < 12 * SUDOKU_GRID_SIZE; i++)
    printf( "%c", (i % 12) ? DOUBLE_FRAME_HORIZONTAL_BAR : DOUBLE_FRAME_TOP_JUNCTION );
  printf( "%c\n", DOUBLE_FRAME_UR_CORNER );

  int i = 0;
  for (int row = 0; row < SUDOKU_GRID_SIZE; row++) {
    printf( "%c", DOUBLE_FRAME_VERTICAL_BAR );
    for (int col = 0; col < SUDOKU_GRID_SIZE; col++, i++) {
      const auto &cell = (*grid_)[i];
      if (cell.solved())
        printf( "   : %d :   %c",
                cell.clue_,
                ((col + 1) % SUDOKU_BOX_SIZE == 0) ? DOUBLE_FRAME_VERTICAL_BAR : SINGLE_FRAME_VERTICAL_BAR );
      else
        printf( " %c%c%c%c%c%c%c%c%c %c",
                (cell.candidate_.clues_[0] != 0) ? '1' : '-',
                (cell.candidate_.clues_[1] != 0) ? '2' : '-',
                (cell.candidate_.clues_[2] != 0) ? '3' : '-',
                (cell.candidate_.clues_[3] != 0) ? '4' : '-',
                (cell.candidate_.clues_[4] != 0) ? '5' : '-',
                (cell.candidate_.clues_[5] != 0) ? '6' : '-',
                (cell.candidate_.clues_[6] != 0) ? '7' : '-',
                (cell.candidate_.clues_[7] != 0) ? '8' : '-',
                (cell.candidate_.clues_[8] != 0) ? '9' : '-',
                ((col + 1) % SUDOKU_BOX_SIZE == 0) ? DOUBLE_FRAME_VERTICAL_BAR : SINGLE_FRAME_VERTICAL_BAR );
    }
    printf( "\n" );
    if (row == SUDOKU_GRID_SIZE - 1)
      continue;
    printf( "%c", DOUBLE_FRAME_LEFT_JUNCTION );
    for (int j = 1; j < 12 * SUDOKU_GRID_SIZE; j++)
      if ((row + 1) % SUDOKU_BOX_SIZE != 0) {
        printf( "%c", (j % 12) ? SINGLE_FRAME_HORIZONTAL_BAR :
                      (j % (12 * SUDOKU_BOX_SIZE)) ? SINGLE_FRAME_CROSS : DOUBLE_FRAME_CROSS );
      } else
        printf( "%c", (j % 12) ? DOUBLE_FRAME_HORIZONTAL_BAR : /*DOUBLE_FRAME_HORIZONTAL_BAR*/ DOUBLE_FRAME_CROSS );
    printf( "%c\n", DOUBLE_FRAME_RIGHT_JUNCTION );
  }

  printf( "%c", DOUBLE_FRAME_LL_CORNER );                  // bottom frame
  for (int i = 1; i < 12 * SUDOKU_GRID_SIZE; i++)
    printf( "%c", (i % 12) ? DOUBLE_FRAME_HORIZONTAL_BAR : DOUBLE_FRAME_BOTTOM_JUNCTION );
  printf( "%c\n", DOUBLE_FRAME_LR_CORNER );

}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::sudokuFinished( void ) {
  // return true, if all clues is sudoku grid are solved

  for (auto &cell : *grid_)
    if (cell.unsolved())
      return false;

  return true;
}  // -------------------------------------------------------------------------------------------

bool SudokuGrid::sudokuSolvingCheck( void ) {
  // return true, if all clues is sudoku grid are correct

  for (auto &cell : *grid_) {
    if (cell.unsolved())
      return false;

    auto currentClue = cell.getClueIndex() + 1;
    cell.resetClue();                                      // because validClue() implementation
    if (validClue( cell, currentClue ) == false) {
      cell.setClue( currentClue );
      return false;
    }
    cell.setClue( currentClue );
  }

  return true;
}  // -------------------------------------------------------------------------------------------

int SudokuGrid::sudokuSolve( int riddleId ) {
  Analysis analyses[] = {
    {&SudokuGrid::singleSolution,   "singleSolution",   /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::singleCell,       "singleCell",       /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::singleBox,        "singleBox",        /* count: */ 0, /* verbose level: */ 0},
  #if 1
    {&SudokuGrid::nakedDuo,         "nakedDuo",         /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedTrio,        "nakedTrio",        /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedQuartet,     "nakedQuartet",     /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedQuintet,     "nakedQuintet",     /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedSextet,      "nakedSextet",      /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedSeptet,      "nakedSeptet",      /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::nakedOctet,       "nakedOctet",       /* count: */ 0, /* verbose level: */ 0},
  #endif
  //{&SudokuGrid::disjoinSubset,    "disjoinSubset",    /* count: */ 0, /* verbose level: */ 0},
  //{&SudokuGrid::disjoinChain,     "disjoinChain",     /* count: */ 0, /* verbose level: */ 0},
  #if 1
    {&SudokuGrid::x_Wing,           "X-Wing",           /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::swordfish,        "swordfish",        /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::jellyfish,        "jellyfish",        /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::starfish,         "starfish",         /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::sixfish,          "sixfish",          /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::sevenfish,        "sevenfish",        /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::octopus,          "octopus",          /* count: */ 0, /* verbose level: */ 0},
  #endif
  #if 1
    {&SudokuGrid::w_Wing,           "W-Wing",           /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::xy_Wing,          "XY-Wing",          /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::xyz_Wing,         "XYZ-Wing",         /* count: */ 0, /* verbose level: */ 0},
  #endif
  #if 1
    {&SudokuGrid::skyscraper,       "skyscraper",       /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::twoStringKite,    "2-String Kite",    /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::emptyRectangle,   "empty rectangle",  /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::finnedX_Wing,     "Finned X-Wing",    /* count: */ 0, /* verbose level: */ 0},
  #endif
  #if 1
    {&SudokuGrid::uniqueRectangle,  "unique rectangle", /* count: */ 0, /* verbose level: */ 0},   // UR type 1...7
  #endif
  //{&SudokuGrid::x_Chain,          "X-Chain",          /* count: */ 0, /* verbose level: */ 0},
    {&SudokuGrid::x_NodeChain,      "X-NodeChain",      /* count: */ 0, /* verbose level: */ 0},
  //{&SudokuGrid::remotePair,       "remote pair",      /* count: */ 0, /* verbose level: */ 0},   // todo!
  #if 0
    {&SudokuGrid::sortedBruteForce, "sortedBruteForce", /* count: */ 0, /* verbose level: */ 0},   // recursive. slow
  #elif 0
    {&SudokuGrid::bruteForce,       "bruteForce",       /* count: */ 0, /* verbose level: */ 0},   // recursive, slowest
  #endif
  };

  printf( "\nsudokuSolve (id: %d, \"%s\")\n", riddleId, riddleName_ );
  initGrid();
  for (bool progressing = true; progressing; ) {
    progressing = false;
    for (auto &analysis : analyses)
      if ((this->*analysis.solvingRecipe)()) {
        printSudokuGrid( "--> %s() successful\n", analysis.name );
        progressing = true;
        break;
      }
  }

  if (sudokuFinished() == false) {
    printSudokuGrid( "\n!!! Solving NOT finished\n" );
    return 0;
  } else if (sudokuSolvingCheck()) {
    printSudokuGrid( "Solving check: Ok\n" );
    return 1;
  } else {
    printSudokuGrid( "\n!!! Solving finished, grid in ERROR !!!\n" );
    return -1;
  }
}  // -------------------------------------------------------------------------------------------

