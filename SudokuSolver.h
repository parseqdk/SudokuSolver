//
// Sudoku Solver - engine
// by d/x, Spring/Summer/Autumn 2022, Daniel Koziarski
// source indentation = 2.
//

//#pragma once
//#pragma pack(1)  // disables padding structures with hidden bytes for 32/64 bit alignment

#ifndef _SUDOKU_SOLVER_H
#define _SUDOKU_SOLVER_H

//#define NDEBUG

#include <cassert>
#include <cstddef>
#include <iterator>
#include <array>
#include <map>
#include <algorithm>
#include <cstring>
//#include <cstdio>
//#include <vector>
//#include <stack>

#if 0
// todo! test template versions od elementsof and capacityof
template<typename T>                                       // C6384 when used with pointer declared like ...(*ptr)[N];
constexpr auto elementsof(T x) { return (sizeof (x) / sizeof ((x)[0])); }

template<typename T>
constexpr auto capacityof(T x) { return static_cast<size_t>( (&(x))[1] - (x) ); }
#else
#define elementsof(x)  (sizeof (x) / sizeof ((x)[0]))      // C6384 when used with pointer declared like ...(*ptr)[N];
#define capacityof(x)  (static_cast<size_t>( (&(x))[1] - (x) ))
#endif

#define last_elementof(x)  ((x) + elementsof(x))

#if 0
constexpr auto SUDOKU_GRID_SIZE = 9    // rows _AND_ columns;
constexpr auto SUDOKU_BOX_SIZE  = 3;
#else
#define SUDOKU_GRID_SIZE  9            // rows _AND_ columns
#define SUDOKU_BOX_SIZE   3
#endif

#define DOUBLE_FRAME_UL_CORNER         201
#define DOUBLE_FRAME_UR_CORNER         187
#define DOUBLE_FRAME_LL_CORNER         200
#define DOUBLE_FRAME_LR_CORNER         188
#define DOUBLE_FRAME_HORIZONTAL_BAR    205
#define DOUBLE_FRAME_VERTICAL_BAR      186
#define DOUBLE_FRAME_LEFT_JUNCTION     204
#define DOUBLE_FRAME_RIGHT_JUNCTION    185
#define DOUBLE_FRAME_TOP_JUNCTION      203
#define DOUBLE_FRAME_BOTTOM_JUNCTION   202
#define DOUBLE_FRAME_CROSS             206

#define SINGLE_FRAME_UL_CORNER         218
#define SINGLE_FRAME_UR_CORNER         191
#define SINGLE_FRAME_LL_CORNER         192
#define SINGLE_FRAME_LR_CORNER         217
#define SINGLE_FRAME_HORIZONTAL_BAR    196
#define SINGLE_FRAME_VERTICAL_BAR      179
#define SINGLE_FRAME_LEFT_JUNCTION     195
#define SINGLE_FRAME_RIGHT_JUNCTION    180
#define SINGLE_FRAME_TOP_JUNCTION      194
#define SINGLE_FRAME_BOTTOM_JUNCTION   193
#define SINGLE_FRAME_CROSS             197

typedef unsigned char uchar;

class SudokuGrid {
  struct Possibilities;
  class Cell {
  protected:
    uchar clue_ {};

    Cell *rowOfCellsStart_ {};
    Cell *colOfCellsStart_ {};
    Cell *boxOfCellsStart_ {};

    struct CandidateTraits {           // todo! replace struct by class ?
      uchar crossOutProtected_ { 0 };
      uchar clues_[SUDOKU_GRID_SIZE] {};

      bool operator == ( const CandidateTraits &rhs ) const {
        for (int i = 0; i < elementsof( this->clues_ ); i++)
          if (this->clues_[i] != rhs.clues_[i])
            return false;
        return true;
      }
      bool operator != ( const CandidateTraits &rhs ) const { return (*this == rhs) == false; }

    #if 1  // 13.08.2022
      friend CandidateTraits operator + ( CandidateTraits lhs, const CandidateTraits &rhs ) {
        for (size_t i = 0; i < elementsof( lhs.clues_ ); i++)
          lhs.clues_[i] = (lhs.clues_[i] != 0 || rhs.clues_[i] != 0) ? 1 : 0;
        return lhs;
      }
      friend CandidateTraits operator - ( CandidateTraits lhs, const CandidateTraits &rhs ) {
        for (size_t i = 0; i < elementsof( lhs.clues_ ); i++)
          lhs.clues_[i] = (lhs.clues_[i] > rhs.clues_[i]) ? lhs.clues_[i] - rhs.clues_[i] : 0;
        return lhs;
      }
      friend CandidateTraits operator * ( CandidateTraits lhs, const CandidateTraits &rhs ) {
        for (size_t i = 0; i < elementsof( lhs.clues_ ); i++)
          lhs.clues_[i] = (lhs.clues_[i] != 0 && rhs.clues_[i] != 0) ? 1 : 0;
        return lhs;
      }
    #else
      CandidateTraits operator + ( const CandidateTraits &rhs ) const {
        CandidateTraits result;
        for (size_t i = 0; i < elementsof( result.clues_ ); i++)
          result.clues_[i] = (this->clues_[i] != 0 || rhs.clues_[i] != 0) ? 1 : 0;
        return result;
      }
      CandidateTraits operator - ( const CandidateTraits &rhs ) const {
        CandidateTraits result;
        for (size_t i = 0; i < elementsof( result.clues_ ); i++)
          result.clues_[i] = (this->clues_[i] > rhs.clues_[i]) ? this->clues_[i] - rhs.clues_[i] : 0;
        return result;
      }
      CandidateTraits operator * ( const CandidateTraits &rhs ) const {
        CandidateTraits result;
        for (size_t i = 0; i < elementsof( result.clues_ ); i++)
          result.clues_[i] = (this->clues_[i] != 0 && rhs.clues_[i] != 0) ? 1 : 0;
        return result;
      }
    #endif

      friend CandidateTraits operator + ( CandidateTraits lhs, int rhsClueIndex ) { return lhs += rhsClueIndex; }
      friend CandidateTraits operator - ( CandidateTraits lhs, int rhsClueIndex ) { return lhs -= rhsClueIndex; }
      CandidateTraits &operator += ( int rhsClueIndex ) { this->clues_[rhsClueIndex] = 1; return *this; }
      CandidateTraits &operator -= ( int rhsClueIndex ) { this->clues_[rhsClueIndex] = 0; return *this; }
      operator int() const { return this->candidateFirstIndex(); }
      CandidateTraits &operator = ( const CandidateTraits &rhs ) {             /* shallow copy is ok */ 
      #if 1
        std::memcpy( this, &rhs, sizeof( *this ) );
      #else
        crossOutProtected_ = rhs.crossOutProtected_;
        *clues_            = *rhs.clues_;
      #endif
        return *this;
      }

      int candidateCount( void ) const {
        int count = 0;
        for (const auto it : clues_)
          count += (it != 0) ? 1 : 0;
        return count;
      }  // -------------------------------------------------------------------------------------------
      bool traitsProtected  ( void ) const { return (crossOutProtected_) ? true : false; }
      bool traitsUnprotected( void ) const { return (crossOutProtected_) ? false : true; }
      bool candidatePossible  ( int clueIndex ) const { return (clues_[clueIndex] != 0) ? true : false; }
      bool candidateImpossible( int clueIndex ) const { return (clues_[clueIndex] == 0) ? true : false; }
      int  candidateNextIndex( int previousIndex ) const {
        // return candidate's next active clueIndex, following previousIndex

        for (int clueIndex = previousIndex + 1; clueIndex < elementsof( clues_ ); clueIndex++)
          if (candidatePossible( clueIndex ))
            return clueIndex;

        return -1;           // caller's algrithm should never lead to this return value
      }  // -------------------------------------------------------------------------------------
      int  candidateFirstIndex( void ) const { return candidateNextIndex( -1 ); }
      //~CandidateTraits() {}
      CandidateTraits() {}
      CandidateTraits( int candidate0, int candidate1 ) {
      #if 1
        crossOutProtected_ = false;

        for (auto &item : clues_)
          item = 0;
      #endif
        clues_[candidate0] = clues_[candidate1] = 1;
      }
    } candidate_;

    int  candidateNextIndex( int previousIndex ) const { return candidate_.candidateNextIndex( previousIndex );}
    int  candidateFirstIndex( void ) const { return candidate_.candidateNextIndex( -1 ); }

    int  candidateCount( void ) const { assert( this->unsolved() ); return this->candidate_.candidateCount(); }

    void candidateProtection ( bool state ) const {
      //const_cast<uchar &> (candidate_.crossOutProtected_) = state;
      assert( state != false && candidate_.crossOutProtected_ < 2  ||
              state == false && candidate_.crossOutProtected_ > 0 );

      const_cast<uchar &> (candidate_.crossOutProtected_) += (state) ? 1 : -1;
    }
    bool candidateProtected  ( void ) const { return (candidate_.crossOutProtected_ != 0); }
    bool candidateUnprotected( void ) const { return (candidate_.crossOutProtected_ == 0); }

    bool candidatePossible  ( int clueIndex ) const { return (candidate_.clues_[clueIndex] != 0); }
    bool candidateImpossible( int clueIndex ) const { return (candidate_.clues_[clueIndex] == 0); }

    void candidateCrossOut( int clueIndex ) { candidate_.clues_[clueIndex] = 0; }
    bool candidateCrossOut( const CandidateTraits &candidateTraits );

    bool candidateCrossOutInProtectedCell( int clueIndex );
    bool candidateCrossOutInProtectedCell( const CandidateTraits &candidateTraits );

    bool candidateFitsInto( const Cell *cell ) const;
    bool candidateMisfitsInto( const Cell *cell ) const { return (candidateFitsInto( cell ) == false); };

    bool candidateIdentical( const Cell *cell ) const { return this->candidate_ == cell->candidate_; }
    bool candidateDifferent( const Cell *cell ) const { return this->candidate_ != cell->candidate_; }

    bool candidateIdentical( const CandidateTraits &candidateTraits ) const { return this->candidate_ == candidateTraits; }
    bool candidateDifferent( const CandidateTraits &candidateTraits ) const { return this->candidate_ != candidateTraits; }

    bool solved  ( void ) const { return (clue_ != 0) ? true : false; }
    bool unsolved( void ) const { return (clue_ == 0) ? true : false; }
    void setClue( int clue ) { clue_ = static_cast<uchar> (clue); }
    void resetClue( void ) { clue_ = 0; }
    int  getClueIndex( void ) const { return (int) clue_ - 1; }

    void mergeCellPossibilities( Possibilities *mergedSum ) const {
      if (unsolved())                  // test for speed
        *mergedSum += *this;           // merge Cell address and possibility infos in Possibilities
    } // ----------------------------------------------------------------------------------------

    bool sharingSameRow( const Cell *cell ) const {
      assert( this != cell );
      return (rowOfCellsStart_ == cell->rowOfCellsStart_);
    } // ----------------------------------------------------------------------------------------
    bool sharingSameCol( const Cell *cell ) const {
      assert( this != cell );
      return (colOfCellsStart_ == cell->colOfCellsStart_);
    } // ----------------------------------------------------------------------------------------
    bool sharingSameBox( const Cell *cell ) const {
      assert( this != cell );
      return (boxOfCellsStart_ == cell->boxOfCellsStart_);
    } // ----------------------------------------------------------------------------------------
    bool sharingSameBox_sameEntityAllowed( const Cell *cell ) const {
      return (boxOfCellsStart_ == cell->boxOfCellsStart_);
    } // ----------------------------------------------------------------------------------------
    bool sharingSameHouse( const Cell *cell ) const {
      assert( this != cell );
      // this and cell are different objects and belong to the same house (row/col/box)
      return (rowOfCellsStart_ == cell->rowOfCellsStart_ ||
              colOfCellsStart_ == cell->colOfCellsStart_ ||
              boxOfCellsStart_ == cell->boxOfCellsStart_);
    } // ----------------------------------------------------------------------------------------
    bool occupyingDifferentRow( const Cell *cell )   const { return (sharingSameRow( cell )   == false); }
    bool occupyingDifferentCol( const Cell *cell )   const { return (sharingSameCol( cell )   == false); }
    bool occupyingDifferentBox( const Cell *cell )   const { return (sharingSameBox( cell )   == false); }
    bool occupyingDifferentHouse( const Cell *cell ) const { return (sharingSameHouse( cell ) == false); }

    bool inAllNeighboursSight( const Cell *const *neighbours, size_t neighboursQuantity ) const {
      // return true, if current cell object is in sight of all neighbours

      for (; neighboursQuantity-- > 0; ) {
        const Cell *neighbour = *neighbours++;
        if (rowOfCellsStart_ != neighbour->rowOfCellsStart_  &&
            colOfCellsStart_ != neighbour->colOfCellsStart_  &&
            boxOfCellsStart_ != neighbour->boxOfCellsStart_)
          return false;                // current object does NOT share the same row/col/box with one of the neighbours
      }
      return true;
    } // ----------------------------------------------------------------------------------------

    int  singleSolution( void ) const;
    void cellInit( uchar clue, int row, int col, std::array<Cell, SUDOKU_GRID_SIZE *SUDOKU_GRID_SIZE> *grid );

  public:
    //~Cell() {}
    //Cell() {}

    friend class SudokuGrid;
  };  // ----------------------------------------------------------------------------------------

  struct Possibilities {
    struct {
      uchar count_;
      Cell *cell_;
    } possibilities_[SUDOKU_GRID_SIZE] = {0};

    Possibilities &operator += ( const Cell &rhs );

    void reset( void ) {
      for (auto &possibility : possibilities_) {
        possibility.count_ = 0;
      //possibility.cell_  = nullptr;  intentionally w/o initialization, cell_ pointer is irrelevant if count_ == 0
      }
    }
  };  // ----------------------------------------------------------------------------------------

  typedef struct {
    bool      (SudokuGrid::*solvingRecipe)( void );
    const char *name;        // for protocol purpose
    size_t     count;        // final statistic
    char       verbose;      // amount of console output protocol
  } Analysis;  // -------------------------------------------------------------------------------

  class HouseIterator;

  template<size_t N>
  class NakedMethods {       // nakedDuo, nakedTriple, nakedQuartet, n..Quintet, n..Sextet, n..Septet, n..Octet
    Cell *naked_[N] {};

    template<size_t N, size_t Step>
    bool analysisStepURtype3( HouseIterator &houseIterator, int houseIndex, const Cell::CandidateTraits &collectedTraits ) {
      bool found = false;
      static_assert( Step <= N - 1  &&  Step > 0 );

      for (; houseIndex < SUDOKU_GRID_SIZE - (Step - 1); houseIndex++) {
        const Cell *cell = naked_[N - Step] = &houseIterator[houseIndex];
        if (cell->solved() || cell->candidateProtected()|| cell->candidateCount() > N)
          continue;
        Cell::CandidateTraits joinedTraits = cell->candidate_ + collectedTraits;
        if (joinedTraits.candidateCount() > N)
          continue;

        if (analysisStepURtype3<N, Step - 1>( houseIterator, houseIndex + 1, joinedTraits )) {
          found = true;
          //break;           /* don't break to proceed with searching, because further finding is indeed possible */
        }
      }
      return found;
    }  // ---------------------------------------------------------------------------------------

    template<size_t N, size_t Step>
    bool analysisStep( HouseIterator &houseIterator, int houseIndex, const Cell::CandidateTraits &collectedTraits ) {
      bool found = false;
      static_assert( Step <= N - 1  &&  Step > 0 );

      for (; houseIndex < SUDOKU_GRID_SIZE - (Step - 1); houseIndex++) {
        const Cell *cell = naked_[N - Step] = &houseIterator[houseIndex];
        if (cell->solved() || cell->candidateCount() > N)
          continue;
        Cell::CandidateTraits joinedTraits = cell->candidate_ + collectedTraits;
        if (joinedTraits.candidateCount() > N)
          continue;

        if (analysisStep<N, Step - 1>( houseIterator, houseIndex + 1, joinedTraits )) {
          found = true;
          //break;           /* don't break to proceed with searching, because further finding is indeed possible */
        }
      }
      return found;
    }  // ---------------------------------------------------------------------------------------

    template<>
    bool analysisStepURtype3<N, 0>( HouseIterator &houseIterator, int houseIndex, const Cell::CandidateTraits &traits ) {
      // This template specialization lets the compiler see, that template's recursion will
      // stop when <Step> reaches a particular value (here: zero).
      // It avoids the C1202 (MS Visual Studio) compiler error: "recursive type or dependency too complex".
      // 
      // (<Step> template argument runs downwards, so it will reach zero.
      //   If <Step> would be running upwards, then a template specialization for a particular <Step>
      //   value greater than start value should also avoid the C1202 error.)
      bool found = false;
      // N-naked constellation found! Showdown: try to cross out.
      chainProtection( true, &naked_[1], N - 1);           // skip unused naked_[0] because of virtual Cell in UR type 3
      if (candidateCrossOutInProtectedHouse( houseIterator, traits ))
        found = true;
      chainProtection( false, &naked_[1], N - 1);

      return found;
    }  // ---------------------------------------------------------------------------------------

    template<>
    bool analysisStep<N, 0>( HouseIterator &houseIterator, int houseIndex, const Cell::CandidateTraits &traits ) {
      // This template specialization lets the compiler see, that template's recursion will
      // stop when <Step> reaches a particular value (here: zero).
      // It avoids the C1202 (MS Visual Studio) compiler error: "recursive type or dependency too complex".
      // 
      // (<Step> template argument runs downwards, so it will reach zero.
      //   If <Step> would be running upwards, then a template specialization for a particular <Step>
      //   value greater than start value should also avoid the C1202 error.)
      bool found = false;
      // N-naked constellation found! Showdown: try to cross out.
      chainProtection( true, naked_, N );
      if (candidateCrossOutInProtectedHouse( houseIterator, traits ))
        found = true;
      chainProtection( false, naked_, N );

      return found;
    }  // ---------------------------------------------------------------------------------------

  #if 1
    bool houseAnalysisURtype3( HouseIterator &&houseIterator, const Cell::CandidateTraits &virtualCellTraits ) {
      // virtualCellTraits are made from third and fourth corner cells of unique rectangle type 3
      // hence, naked_[0] location remains unused
      bool found = false;

      if (virtualCellTraits.candidateCount() <= N)
        if (analysisStepURtype3<N, N - 1>( houseIterator, 0, virtualCellTraits ))
          found = true;

      return found;
    }  // ---------------------------------------------------------------------------------------

    bool houseAnalysis( HouseIterator &houseIterator ) {
      bool found = false;

      for (int houseIndex = 0; houseIndex < SUDOKU_GRID_SIZE - (N - 1); houseIndex++) {
        const Cell *cell = naked_[0] = &houseIterator[houseIndex];
        if (cell->solved() || cell->candidateCount() > N)
          continue;

        if (analysisStep<N, N - 1>( houseIterator, houseIndex + 1, cell->candidate_ )) {
          found = true;
          //break;           /* don't break to proceed with searching, because further finding is indeed possible */
        }
      }
      return found;
    }  // ---------------------------------------------------------------------------------------
  public:
    bool analysisURtype3( Cell *third, Cell *fourth, const Cell::CandidateTraits &virtualCellTraits ) {
      // return false, if neither clue possibility can be crossed out
      // search as usual for naked constellation along houses shared by third and fourth, but:
      //   - ignore third and fourth unique rectangle corner cells (which are made protected)
      //   - use virtualCellTraits made from third and fourth corner cell candidates
      bool found = false;

      third->candidateProtection( true ), fourth->candidateProtection( true );

      if (third->sharingSameBox( fourth ))
        if (houseAnalysisURtype3( BoxIterator { third->boxOfCellsStart_ }, virtualCellTraits ))
          found = true;
      if (third->sharingSameRow( fourth ))
        if (houseAnalysisURtype3( RowIterator { third->rowOfCellsStart_ }, virtualCellTraits ))
          found = true;
      if (third->sharingSameCol( fourth ))
        if (houseAnalysisURtype3( ColIterator { third->colOfCellsStart_ }, virtualCellTraits ))
          found = true;

      third->candidateProtection( false ), fourth->candidateProtection( false );

      return found;
    }  // ---------------------------------------------------------------------------------------

    bool analysis( HouseIterator &&houseIterator ) {
      // return false, if neither clue possibility can be crossed out
      // N == 2 : naked duo (aka: naked pair) algorithm level: intermediate. Description:
      //   YouTube: "dxSudoku #4 Naked Pair"
      //
      //   Obviously possible type of naked duo: XY, XY.
      // 
      // N == 3 : naked triple algorithm level: intermediate. Description:
      //   YouTube: "dxSudoku #7, #8 Naked Triple"
      // 
      //   Possible type of naked triple: XYZ-XYZ-XYZ, XYZ-XYZ-XY, XYZ-XZ-YZ, XY-XZ-YZ (all types in any combination)
      // 
      // N == 4 : naked quartet (aka: naked quad) algorithm level: advanced. Description:
      //   YouTube: "dxSudoku #57 Naked Quad"
      // 
      //   Possible type of naked quartet: WXYZ-WXYZ-WXYZ-WXYZ, ....., WX-XY-YZ-WZ (all types in any combination)
      // 
      // N == 5 : naked quintet algorithm level: extreme. Description: 1 step  further than naked quartet.
      //   Youtube: "dxSudoku #81 Hidden Quadruples and Naked Quintets"
      // N == 6 : naked sextet  algorithm level: ?. Description: 2 steps further than naked quartet.
      // N == 7 : naked septet  algorithm level: ?. Description: 3 steps further than naked quartet.
      // N == 8 : naked octet   algorithm level: ?. Description: 4 steps further than naked quartet.
      // 
      //   Same algorithm should be iterated row by row, then column by column, and then box by box.
      bool found = false;

      for (int house = 0; house < SUDOKU_GRID_SIZE; house++, houseIterator.nextHouse())
        if (houseAnalysis( houseIterator ))
          found = true;

      return found;
    }  // ---------------------------------------------------------------------------------------
  #else
  public:
    bool analysis( HouseIterator &&houseIterator ) {
      // return false, if neither clue possibility can be crossed out
      // N == 2 : naked duo (aka: naked pair) algorithm level: intermediate. Description:
      //   YouTube: "dxSudoku #4 Naked Pair"
      //
      //   Obviously possible type of naked duo: XY, XY.
      // 
      // N == 3 : naked triple algorithm level: intermediate. Description:
      //   YouTube: "dxSudoku #7, #8 Naked Triple"
      // 
      //   Possible type of naked duo: XYZ-XYZ-XYZ, XYZ-XYZ-XY, XYZ-XZ-YZ, XY-XZ-YZ (all types in any combination)
      // 
      // N == 4 : naked quartet (aka: naked quad) algorithm level: advanced. Description:
      //   YouTube: "dxSudoku #57 Naked Quad"
      // 
      //   Possible type of naked duo: WXYZ-WXYZ-WXYZ-WXYZ, ....., WX-XY-YZ-WZ (all types in any combination)
      // 
      // N == 5 : naked quintet algorithm level: extreme. Description: 1 step  further than naked quartet.
      //   Youtube: "dxSudoku #81 Hidden Quadruples and Naked Quintets"
      // N == 6 : naked sextet  algorithm level: ?. Description: 2 steps further than naked quartet.
      // N == 7 : naked septet  algorithm level: ?. Description: 3 steps further than naked quartet.
      // N == 8 : naked octet   algorithm level: ?. Description: 4 steps further than naked quartet.
      // 
      //   Same algorithm should be iterated row by row, then column by column, and then box by box.
      bool found = false;

      for (int house = 0; house < SUDOKU_GRID_SIZE; house++, houseIterator.nextHouse())
        for (int cellIndex = 0; cellIndex < SUDOKU_GRID_SIZE - (N - 1); cellIndex++) {
          const Cell *cell = naked_[0] = &houseIterator[cellIndex];
          if (cell->solved() || cell->candidateCount() > N)
            continue;

          if (analysisStep<N, N - 1>( houseIterator, cellIndex + 1, cell->candidate_ )) {
            found = true;
            //break;         /* don't break to carry on searching, because further finding is indeed possible */
          }
        }
      return found;
    }  // ---------------------------------------------------------------------------------------
  #endif
    friend class SudokuGrid;
  };
  using NakedMethodDuo     = NakedMethods<2>;
  using NakedMethodTrio    = NakedMethods<3>;
  using NakedMethodQuartet = NakedMethods<4>;
  using NakedMethodQuintet = NakedMethods<5>;
  using NakedMethodSextet  = NakedMethods<6>;
  using NakedMethodSeptet  = NakedMethods<7>;
  using NakedMethodOctet   = NakedMethods<8>;

  template<size_t N> class FishMethods;
  template<size_t N>
  class CellFishSequence {   // x-Wing, swordfish, jellyfish, starfish, sixfish, sevenfish, octopus
    bool rowMode_;
    Cell *fish_[N];

    CellFishSequence() : rowMode_ { false }, fish_ { nullptr } {}
    CellFishSequence( HouseIterator &houseIterator, bool rowMode, int clueIndex );
    CellFishSequence( HouseIterator &houseIterator, bool rowMode, Possibilities &possibilities, int clueIndex );
    // todo! dlaczego jest constructor CellFishSequence() i CellFishSequenceInit() o identycznej zawartoœci?
    // todo! jeœli faktycznie potrzebne s¹ oba warianty, to constructor powininen wywo³ywaæ CellFishSequenceInit()
    //       i nie powtarzaæ tego samego kodu!
    void CellFishSequenceInit( bool rowMode, Possibilities &possibilities, int clueIndex);

    static bool invalidClueCount( uchar clueCount ) { return ( clueCount < 2 || clueCount > N ); }

    void chainProtection( bool state ) {
      for (size_t i = 0; i < elementsof( fish_ ); i++) {
        assert( fish_[i] != nullptr );
        fish_[i]->candidateProtection( state );
      }
    }  // ---------------------------------------------------------------------------------------

    static void chainProtection( bool state, CellFishSequence<N> (&sequence)[N], size_t nSequences = N ) {
      // update candidate protection state for all cells in sequences (N * N)
      for (size_t n = 0; n < nSequences; n++)
        for (size_t i = 0; i < elementsof( fish_ ); i++)
          if (sequence[n].fish_[i] != nullptr)
            sequence[n].fish_[i]->candidateProtection( state );
    }  // ---------------------------------------------------------------------------------------

    size_t occupiedBoxesCount( void ) {
      // return: quantity of different sudoku's boxes occupied by the elements in this->fish_[N];
      // !! complexity: N*(N-1)/2 (worst case), i.e. acceptable for moderate N
      size_t uniqueBoxCount = 0;
    #if 0          // correct version (STL:std::none_of() in std::for_each() with nested lambdas)
      auto last = fish_;
      auto freshBoxCounter = [&]( const auto cell ) {
        const auto freshBox = cell->boxOfCellsStart_;
        uniqueBoxCount += std::none_of( fish_, last++, [&]( const auto it ) { return it->boxOfCellsStart_ == freshBox; } );
      };
      std::for_each( fish_, fish_ + elementsof( fish_ ), freshBoxCounter );
    #elif 1        // correct version (STL:std::none_of() in range loop)
      auto last = fish_;
      for (const auto cell : fish_) {
        const auto freshBox = cell->boxOfCellsStart_;
        uniqueBoxCount += std::none_of( fish_, last++, [&]( const auto it ) { return it->boxOfCellsStart_ == freshBox; } );
      }
    #elif 1        // correct version (STL:std::none_of() in counting loop )
      for (size_t i = 0; i < elementsof( fish_ ); i++) {
        const auto freshBox          = fish_[i]->boxOfCellsStart_;
        const auto sameBoxComparator = [&]( const auto it ) { return it->boxOfCellsStart_ == freshBox; };
        uniqueBoxCount += std::none_of( fish_, fish_ + i, sameBoxComparator );
      }
    #else          // correct version (without STL)
      auto uniqueBox = []( Cell **first, Cell **const last) -> size_t {
        const auto freshBox = (*last)->boxOfCellsStart_;
        while (first != last)
          if ((*first++)->boxOfCellsStart_ == freshBox)
            return 0;
        return 1;            // box is NOT in the range [*first, last), i.e. box is unique
      };
      for (size_t i = 0; i < elementsof( fish_ ); i++)
        uniqueBoxCount += uniqueBox( fish_, fish_ + i );
    #endif
      return uniqueBoxCount;
    }  // ---------------------------------------------------------------------------------------

    bool coversGroup( HouseIterator &groupIterator, int clueIndex ) const {
      // return true, if all group cells (with clueIndex as candidate) are covered by lines of fish_ (i.e. no fin is required)
      for (int groupIndex = 0; groupIndex < SUDOKU_BOX_SIZE; groupIndex++) {
        const Cell &groupCell = groupIterator[groupIndex];
        if (groupCell.solved() || groupCell.candidateImpossible( clueIndex ))
          continue;

        const Cell *groupCellLine = (rowMode_) ? groupCell.colOfCellsStart_ : groupCell.rowOfCellsStart_;
        for (int i = 0; i < elementsof( fish_ ); i++) {
          const Cell *fishCell = fish_[i];
          if (fishCell != nullptr  &&
              groupCellLine == ((rowMode_) ? fishCell->colOfCellsStart_ : fishCell->rowOfCellsStart_))
            break;                                                             // matching found
          else if (i == elementsof( fish_ ) - 1)
            return false;                                                      // mismatch found
        }
      }
      return true;
    }  // ---------------------------------------------------------------------------------------

    bool absentInChute( const Cell *chuteCell ) const {
      // return true, if none of fish_ elements occupies same chute as chuteCell (perpendicular to fish_)

      const Cell *chuteBox = chuteCell->boxOfCellsStart_;
      const Cell *chute = (rowMode_) ? chuteBox->colOfCellsStart_ : chuteBox->rowOfCellsStart_;

      for (const auto cell : fish_) {
        if (cell == nullptr)
          continue;
        const Cell *cellBox = cell->boxOfCellsStart_;
        if (chute == ((rowMode_) ? cellBox->colOfCellsStart_ : cellBox->rowOfCellsStart_))
          return false;
      }
      return true;
    }  // ---------------------------------------------------------------------------------------

    bool presentInChute( const Cell *chuteCell ) const {
      // return true, if any fish_ element occupies same chute as chuteCell (perpendicular to fish_)
      return absentInChute( chuteCell ) == false;
    }  // ---------------------------------------------------------------------------------------

    bool operator != ( const CellFishSequence<N> &rhs ) const {
      Cell  *collector[N] {};
      size_t collectorSize = 0;

      for (size_t i = 0; i < elementsof( collector ); i++)
        collectorSize += ((collector[i] = fish_[i]) != nullptr) ? 1 : 0;

      for (const auto cell : rhs.fish_) {
        if (cell == nullptr)
          continue;

        for (size_t i = 0; i < elementsof( collector ); i++) {
          if (i == collectorSize) {
            collector[collectorSize++] = cell;                                 // append new cell and
            break;                                                             // treat it as matching
          } else if (rowMode_ == true  && cell->colOfCellsStart_ == collector[i]->colOfCellsStart_  ||
                     rowMode_ == false && cell->rowOfCellsStart_ == collector[i]->rowOfCellsStart_)
            break;                                                             // matching found
          else if (i == elementsof( collector ) - 1)
            return true;                                                       // args different
        }
      }
      return false;
    }  // ---------------------------------------------------------------------------------------

    CellFishSequence<N> operator + ( const CellFishSequence<N> &rhs ) const {
      CellFishSequence<N> result = *this;
      size_t                resultSize = 0;

      for (size_t i = 0; i < elementsof( result.fish_ ); i++)
        resultSize += (result.fish_[i] != nullptr) ? 1 : 0;

      for (const auto cell : rhs.fish_) {
        if (cell == nullptr)
          continue;

        for (size_t i = 0; i < elementsof( result.fish_ ); i++) {
          if (i == resultSize) {
            result.fish_[resultSize++] = cell;                                 // append new cell and
            break;                                                             // treat it as matching
          } else if (rowMode_ == true  && cell->colOfCellsStart_ == result.fish_[i]->colOfCellsStart_  ||
                     rowMode_ == false && cell->rowOfCellsStart_ == result.fish_[i]->rowOfCellsStart_)
            break;                                                             // matching found
        }
      }
      return result;
    }  // ---------------------------------------------------------------------------------------

    CellFishSequence<N> &operator += ( const CellFishSequence<N> &rhs ) {
      size_t resultSize = 0;

      for (size_t i = 0; i < elementsof( this->fish_ ); i++)
        resultSize += (this->fish_[i] != nullptr) ? 1 : 0;

      for (const auto cell : rhs.fish_) {
        if (cell == nullptr)
          continue;

        for (size_t i = 0; i < elementsof( this->fish_ ); i++) {
          if (i == resultSize) {
            this->fish_[resultSize++] = cell;                                  // append new cell and
            break;                                                             // treat it as matching
          } else if (rowMode_ == true  && cell->colOfCellsStart_ == this->fish_[i]->colOfCellsStart_  ||
                     rowMode_ == false && cell->rowOfCellsStart_ == this->fish_[i]->rowOfCellsStart_)
            break;                                                             // matching found
        }
      }
      return *this;
    }  // ---------------------------------------------------------------------------------------

    friend class FishMethods<N>;
    friend class SudokuGrid;
  };  // ----------------------------------------------------------------------------------------

  template<> class CellFishSequence<2> {  // specialization for x-Wing
    Cell *fish_[2];
    bool rowMode_;

    CellFishSequence() : fish_ { nullptr }, rowMode_ { false } {}
    CellFishSequence( HouseIterator &houseIterator, bool rowMode, int clueIndex ) : rowMode_ { rowMode }  {
      Cell **fish = fish_;
      size_t fishCapacity = 2;

      for (int i = 0; i < SUDOKU_GRID_SIZE; i++)
        if (houseIterator[i].candidatePossible( clueIndex ) && fishCapacity-- > 0) 
          *fish++ = &houseIterator[i];
    }  // ---------------------------------------------------------------------------------------
    CellFishSequence( HouseIterator &houseIterator, bool rowMode,
                      Possibilities &possibilities, int clueIndex ) : rowMode_ { rowMode } {

      assert( possibilities.possibilities_[clueIndex].count_ == 2 );

      fish_[0] = possibilities.possibilities_[clueIndex].cell_;
      fish_[1] = lookBackForCellWithSamePotentialClue( houseIterator, fish_[0], clueIndex );
    }  // ---------------------------------------------------------------------------------------
    // todo! dlaczego jest constructor CellFishSequence() i CellFishSequenceInit() o identycznej zawartoœci?
    // todo! jeœli faktycznie potrzebne s¹ oba warianty, to constructor powininen wywo³ywaæ CellFishSequenceInit()
    //       i nie powtarzaæ tego samego kodu!
    void CellFishSequenceInit( bool rowMode, Possibilities &possibilities, int clueIndex) {

      rowMode_ = rowMode;

      assert( possibilities.possibilities_[clueIndex].count_ == 2 );

      // make sure possibilities.possibilities_[clueIndex].count_ has been not invalidated (after candidate crossing out):
      assert( possibilities.possibilities_[clueIndex].count_ ==
              candidateCountInHouse( (rowMode)
              ? (HouseIterator &) RowIterator { possibilities.possibilities_[clueIndex].cell_->rowOfCellsStart_ }
              : (HouseIterator &) ColIterator { possibilities.possibilities_[clueIndex].cell_->colOfCellsStart_ }, clueIndex ) );

      fish_[0] = possibilities.possibilities_[clueIndex].cell_;
      HouseIterator &lineIterator = (rowMode) ? (HouseIterator &) RowIterator { fish_[0] } 
                                              : (HouseIterator &) ColIterator { fish_[0] };

      fish_[1] = lookBackForCellWithSamePotentialClue( lineIterator, clueIndex );
    }  // ---------------------------------------------------------------------------------------

    static bool invalidClueCount( uchar clueCount ) { return ( clueCount != 2 ); }

    void chainProtection( bool state ) {
      fish_[0]->candidateProtection( state );
      fish_[1]->candidateProtection( state );
    }  // ---------------------------------------------------------------------------------------

    static void chainProtection( bool state, CellFishSequence<2> (&sequence)[2], size_t nSequences = 2 ) {
      assert( nSequences > 0  &&  nSequences <= 2 );
      for (size_t n = 0; n < nSequences; n++)
        for (size_t i = 0; i < elementsof( fish_ ); i++)
          if (sequence[n].fish_[i] != nullptr)
            sequence[n].fish_[i]->candidateProtection( state );
    }  // ---------------------------------------------------------------------------------------

    bool coversGroup( HouseIterator &groupIterator, int clueIndex ) const {
      // return true, if all group cells (with clueIndex as candidate) are covered by lines of fish_ (i.e. no fin is required)
      for (int groupIndex = 0; groupIndex < SUDOKU_BOX_SIZE; groupIndex++) {
        const Cell &groupCell = groupIterator[groupIndex];
        if (groupCell.solved() || groupCell.candidateImpossible( clueIndex ))
          continue;

        const Cell *groupCellLine = (rowMode_) ? groupCell.colOfCellsStart_ : groupCell.rowOfCellsStart_;
        for (int i = 0; i < elementsof( fish_ ); i++) {
          const Cell *fishCell = fish_[i];
          if (fishCell != nullptr  &&
               groupCellLine == ((rowMode_) ? fishCell->colOfCellsStart_ : fishCell->rowOfCellsStart_))
            break;                                                             // matching found
          else if (i == elementsof( fish_ ) - 1)
            return false;                                                      // mismatch found
        }
      }
      return true;
    }  // ---------------------------------------------------------------------------------------

    bool absentInChute( const Cell *chuteCell ) const {
      // return true, if none of fish_ elements occupies same chute as chuteCell (perpendicular to fish_)

      const Cell *chuteBox = chuteCell->boxOfCellsStart_;
      const Cell *chute = (rowMode_) ? chuteBox->colOfCellsStart_ : chuteBox->rowOfCellsStart_;

      for (const auto cell : fish_) {
        if (cell == nullptr)
          continue;
        const Cell *cellBox = cell->boxOfCellsStart_;
        if (chute == ((rowMode_) ? cellBox->colOfCellsStart_ : cellBox->rowOfCellsStart_))
          return false;
      }
      return true;
    }  // ---------------------------------------------------------------------------------------

    bool presentInChute( const Cell *chuteCell ) const {
      // return true, if any fish_ element occupies same chute as chuteCell (perpendicular to fish_)
      return absentInChute( chuteCell ) == false;
    }  // ---------------------------------------------------------------------------------------

    bool operator != ( const CellFishSequence<2> &rhs ) const {
      for (const auto cell : rhs.fish_) {
        for (size_t i = 0; i < elementsof( fish_ ); i++) {
          if (rowMode_ == true  && cell->colOfCellsStart_ == fish_[i]->colOfCellsStart_  ||
              rowMode_ == false && cell->rowOfCellsStart_ == fish_[i]->rowOfCellsStart_)
            break;                                                             // matching found
          else if (i == elementsof( fish_ ) - 1)
            return true;                                                       // args different
        }
      }
      return false;
    }  // ---------------------------------------------------------------------------------------

    CellFishSequence<2> operator + ( const CellFishSequence<2> &rhs ) const { return *this; }

    CellFishSequence<2> &operator += ( const CellFishSequence<2> &rhs ) { return *this; }

    friend class FishMethods<2>;
    friend class SudokuGrid;
  };  // ----------------------------------------------------------------------------------------

  using CellFishDuo     = CellFishSequence<2>;
  using CellFishTrio    = CellFishSequence<3>;
  using CellFishQuartet = CellFishSequence<4>;
  using CellFishQuintet = CellFishSequence<5>;
  using CellFishSextet  = CellFishSequence<6>;
  using CellFishSeptet  = CellFishSequence<7>;
  using CellFishOctet   = CellFishSequence<8>;

  template<size_t N>
  class FishMethods {
    int                 fishRow_[N] {};
    Possibilities       mergedSum_[SUDOKU_GRID_SIZE];
    CellFishSequence<N> killZone_, fishSequence_[N];

    static bool fishCandidateCrossOut( CellFishSequence<N> (&sequence)[N], const CellFishSequence<N> &killZone, int clueIndex) {
      // return false, if neither clue possibility can be crossed out
      bool found = false;

      CellFishSequence<N>::chainProtection( true, sequence, N );
      for ( const auto &killer : killZone.fish_ )
        if (killer != nullptr) {
          auto &lineIterator = (killZone.rowMode_) ? (HouseIterator &) ColIterator { killer->colOfCellsStart_ }
                                                   : (HouseIterator &) RowIterator { killer->rowOfCellsStart_ };
          if (candidateCrossOutInProtectedHouse( lineIterator, clueIndex ))
            found = true;
        }
      CellFishSequence<N>::chainProtection( false, sequence, N );

      if (found) {
        printf( "%s(clue=%d) in %s: ",
                (N == 2) ? "x-Wing"    :
                (N == 3) ? "swordfish" :
                (N == 4) ? "jellyfish" :
                (N == 5) ? "starfish"  :
                (N == 6) ? "sixfish"   :
                (N == 7) ? "sevenfish" :
                (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (killZone.rowMode_) ? "rows" : "cols" );
        for ( const auto &item : sequence ) {
          for (size_t i = 0; i < N; i++)
            printf( "%s%2d%c", (i == 0) ? "[" : "",
                    (item.fish_[i] == nullptr) ? -1 :
                    (killZone.rowMode_) ? colNumber( item.fish_[i] ) : rowNumber( item.fish_[i] ), (i < N-1) ? ',' : ']' );
        }
        for (size_t i = 0; i < N; i++)
          printf( "%s%2d%s", (i > 0) ? "" : (killZone.rowMode_) ? " kill zone [cols:" : " kill zone [rows:",
                  (killZone.fish_[i] == nullptr) ? -1 :
                  (killZone.rowMode_) ? colNumber( killZone.fish_[i] ) : rowNumber( killZone.fish_[i] ),
                  (i < N-1) ? "," : "] (crossed out: YES)\n" );
      }
      return found;
    }  // ---------------------------------------------------------------------------------------

    static bool finnedFishCandidateCrossOut( Cell *finGroup, CellFishSequence<N> (&sequence)[N],
                                             const CellFishSequence<N> &killZone, int clueIndex ) {
      // return false, if neither clue possibility can be crossed out
      //
      // Candidate eliminations occur in same box as fins if also crossed by lines defined by killZone elements and
      // perpendicular to the killZone sequence.
      // Finned cell group itself and base sets must be excluded from eliminations.
      bool found = false;
      BoxIterator boxIterator { finGroup->boxOfCellsStart_ };
      HouseIterator &finGroupIterator = (killZone.rowMode_) ? (HouseIterator &) RowIterator { finGroup }
                                                            : (HouseIterator &) ColIterator { finGroup };

      CellFishSequence<N>::chainProtection( true, sequence, N - 1 );
      candidateGroupProtection( true, finGroupIterator );
      for ( const auto &killer : killZone.fish_ )
        if (killer != nullptr) {
          for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
            Cell &cell = boxIterator[i];
            if (cell.solved())
              continue;
            if (killZone.rowMode_ == true  && cell.colOfCellsStart_ == killer->colOfCellsStart_  ||
                killZone.rowMode_ == false && cell.rowOfCellsStart_ == killer->rowOfCellsStart_)
              if (cell.candidateCrossOutInProtectedCell( clueIndex ))
                found = true;
          }
        }
      candidateGroupProtection( false, finGroupIterator );
      CellFishSequence<N>::chainProtection( false, sequence, N - 1 );

      return found;
    }  // ---------------------------------------------------------------------------------------

    static bool finnedFishCandidateCrossOut( Cell *finGroup0, Cell *finGroup1, CellFishSequence<N>(&sequence)[N],
                                             const CellFishSequence<N> &killZone, int clueIndex ) {
      // return false, if neither clue possibility can be crossed out
      //
      // Candidate eliminations occur in same box as fins if also crossed by lines defined by killZone elements and
      // perpendicular to the killZone sequence.
      // Finned cell group itself and base sets must be excluded from eliminations.
      bool found = false;
      BoxIterator boxIterator { finGroup0->boxOfCellsStart_ };
      //printf( "finnedFishCandidateCrossOut() clue:%d, TWO fins: @[row:%d, cols:%d],[row:%d, cols:%d]\n",
      //        clueIndex + 1, rowNumber( finGroup0 ), colNumber( finGroup0 ), rowNumber( finGroup1 ), colNumber( finGroup1 ) );
      HouseIterator &finGroupIt0 = (killZone.rowMode_) ? (HouseIterator &) RowIterator { finGroup0 }
                                                       : (HouseIterator &) ColIterator { finGroup0 };
      HouseIterator &finGroupIt1 = (killZone.rowMode_) ? (HouseIterator &) RowIterator { finGroup1 }
                                                       : (HouseIterator &) ColIterator { finGroup1 };
      CellFishSequence<N>::chainProtection( true, sequence, N - 2 );
      candidateGroupProtection( true, finGroupIt0 );
      candidateGroupProtection( true, finGroupIt1 );
      for ( const auto &killer : killZone.fish_ )
        if (killer != nullptr) {
          for (int i = 0; i < SUDOKU_GRID_SIZE; i++) {
            Cell &cell = boxIterator[i];
            if (cell.solved())
              continue;
            if (killZone.rowMode_ == true  && cell.colOfCellsStart_ == killer->colOfCellsStart_  ||
                killZone.rowMode_ == false && cell.rowOfCellsStart_ == killer->rowOfCellsStart_)
              if (cell.candidateCrossOutInProtectedCell( clueIndex ))
                found = true;
          }
        }
      candidateGroupProtection( false, finGroupIt1 );
      candidateGroupProtection( false, finGroupIt0 );
      CellFishSequence<N>::chainProtection( false, sequence, N - 2 );

      return found;
    }  // ---------------------------------------------------------------------------------------

    Cell *finORsashimiLine( HouseIterator &lineIterator, const CellFishSequence<N> &killZone, int clueIndex ) {
      int usualMatch = 0, finMatch = 0;
      Cell *finGroupStart {};

      for (int box = 0; box < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; box++) {
        Cell *groupStart = &lineIterator[box * SUDOKU_BOX_SIZE];
        auto &groupIterator = (killZone.rowMode_) ? (HouseIterator &) RowIterator { groupStart }
                                                  : (HouseIterator &) ColIterator { groupStart };
        if (candidateCountInGroup( groupIterator, clueIndex ) == 0)
          continue;
        if (killZone_.coversGroup( groupIterator, clueIndex ))
          usualMatch += 1;                                 // further fin/sashimi group still possible
        else if (killZone_.absentInChute( &groupIterator[0])) {
          usualMatch = 0;                                  // current line incapable of finned/sashimi N-fish
          break;
        } else if (finMatch != 0) {
          usualMatch = 0;                                  // exactly 1 fin/sashimi per line allowed
          break;
        } else                                             // matching (under fin/sashimi conditions)
          finMatch = 1, finGroupStart = groupStart;
      }
      //if (usualMatch > 0)
      //  printf( "Found fin/sashimi line (%s, clue:%d), usualMatch:%d\n",
      //          (killZone.rowMode_) ? "rows" : "cols", clueIndex + 1, usualMatch );
      return (usualMatch > 0  &&  finMatch == 1) ? finGroupStart : nullptr;
    }  // ---------------------------------------------------------------------------------------

    template<size_t N, size_t Step>
    bool completeFishSequenceWithTwoFins( bool rowMode, int clueIndex ) {
      bool found = false;
    #if 0
      printf( "Chance for TWIN-finned fish: %s clue=%d %s: ",
              (N == 2) ? "x-Wing"    :
              (N == 3) ? "swordfish" :
              (N == 4) ? "jellyfish" :
              (N == 5) ? "starfish"  :
              (N == 6) ? "sixfish"   :
              (N == 7) ? "sevenfish" :
              (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
      for (size_t i = 0; i < N - Step; i++)
        printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
      printf( " (sizes:" );
      for (size_t i = 0; i < N - Step; i++)
        printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
      printf( ") -- _LAST_BUT_TWO_ fishAnalysisStep<%zd, %zu>()\n", N, Step );
    #endif
      Cell *finGroup[SUDOKU_BOX_SIZE] { nullptr };
      Cell *grid = killZone_.fish_[0]->rowOfCellsStart_->colOfCellsStart_;
      HouseIterator &lineIt = (rowMode) ? (HouseIterator &) RowIterator { grid }
                                        : (HouseIterator &) ColIterator { grid };
      for (auto [line, baseSetNo] = std::pair { 0, 0 }; line < SUDOKU_GRID_SIZE; line++, lineIt.nextHouse()) {
        if (baseSetNo < N - Step  &&  line == fishRow_[baseSetNo]) {
          baseSetNo += 1;
          continue;                                        // skip used grid lines during rescanning
        }
        if (&lineIt[0] == lineIt[0].boxOfCellsStart_)      // new chute start?
          for (auto item : finGroup)
            item = nullptr;

        if (finGroup[0] == nullptr)
          finGroup[0] = finORsashimiLine( lineIt, killZone_, clueIndex );
        else if (finGroup[1] == nullptr) {
          finGroup[1] = finORsashimiLine( lineIt, killZone_, clueIndex );
          if (finGroup[1] != nullptr)
            if (finGroup[1]->sharingSameBox( finGroup[0] )) {
            #if 0
              printf( "%s with TWO finned/Sashimi groups, clue:%d %s: ",
                      (N == 2) ? "x-Wing"    :
                      (N == 3) ? "swordfish" :
                      (N == 4) ? "jellyfish" :
                      (N == 5) ? "starfish"  :
                      (N == 6) ? "sixfish"   :
                      (N == 7) ? "sevenfish" :
                      (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
              for (size_t i = 0; i < N - Step; i++)
                printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
              printf( " (sizes:" );
              for (size_t i = 0; i < N - Step; i++)
                printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
               printf( ")\n  Fin/Sashimi found, finGroups @[row:%d, col:%d],[row:%d, col:%d]\n",
                       rowNumber( finGroup[0] ), colNumber( finGroup[0] ), rowNumber( finGroup[1] ), colNumber( finGroup[1] ) );
            #endif
              if (finnedFishCandidateCrossOut( finGroup[0], finGroup[1], fishSequence_, killZone_, clueIndex ))
                found = true;
            }
        } else {
          finGroup[2] = finORsashimiLine( lineIt, killZone_, clueIndex );
          if (finGroup[2] != nullptr) {
            if (finGroup[2]->sharingSameBox( finGroup[0] ))
              if (finnedFishCandidateCrossOut( finGroup[0], finGroup[2], fishSequence_, killZone_, clueIndex))
                found = true;
            else if (finGroup[2]->sharingSameBox( finGroup[1] ))
              if (finnedFishCandidateCrossOut( finGroup[1], finGroup[2], fishSequence_, killZone_, clueIndex))
                found = true;
          }
        }
      }
      return found;
    }  // ---------------------------------------------------------------------------------------

    template<size_t N, size_t Step>
    bool completeFishSequenceWithOneFin( bool rowMode, int clueIndex ) {
      bool found = false;
    #if 0
      printf( "Chance for finned fish: %s clue=%d %s: ",
              (N == 2) ? "x-Wing"    :
              (N == 3) ? "swordfish" :
              (N == 4) ? "jellyfish" :
              (N == 5) ? "starfish"  :
              (N == 6) ? "sixfish"   :
              (N == 7) ? "sevenfish" :
              (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
      for (size_t i = 0; i < N - Step; i++)
        printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
      printf( " (sizes:" );
      for (size_t i = 0; i < N - Step; i++)
        printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
      printf( ") -- _LAST_BUT_ONE_ fishAnalysisStep<%zd, %zu>()\n", N, Step );
    #endif
      Cell *grid = killZone_.fish_[0]->rowOfCellsStart_->colOfCellsStart_;
      HouseIterator &lineIterator = (rowMode) ? (HouseIterator &) RowIterator { grid }
                                              : (HouseIterator &) ColIterator { grid };
      for (auto [line, baseSetNo] = std::pair { 0, 0 }; line < SUDOKU_GRID_SIZE; line++, lineIterator.nextHouse()) {
        if (baseSetNo < N - Step  &&  line == fishRow_[baseSetNo]) {
          baseSetNo += 1;
          continue;                                        // skip used grid lines during rescanning
        }
      #if 1
        Cell *finGroupStart = finORsashimiLine( lineIterator, killZone_, clueIndex );
        if (finGroupStart != nullptr) {
          //HouseIterator &finGroupIterator = (rowMode) ? (HouseIterator &) RowIterator { finGroupStart }
          //                                            : (HouseIterator &) ColIterator { finGroupStart };
          if (finnedFishCandidateCrossOut( finGroupStart, fishSequence_, killZone_, clueIndex )) {
            printf( "%s with ONE finned/Sashimi group, clue:%d %s: ",
                    (N == 2) ? "x-Wing"    :
                    (N == 3) ? "swordfish" :
                    (N == 4) ? "jellyfish" :
                    (N == 5) ? "starfish"  :
                    (N == 6) ? "sixfish"   :
                    (N == 7) ? "sevenfish" :
                    (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
            for (size_t i = 0; i < N - Step; i++)
              printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
            printf( " (sizes:" );
            for (size_t i = 0; i < N - Step; i++)
              printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
            printf( ")\n  Fin/Sashimi found in line:%d, finGroupStart @[row:%d, col:%d]\n",
                    line, rowNumber( finGroupStart ), colNumber( finGroupStart ));
            found = true;
          }
        }
      #else
        int usualMatch = 0, finMatch = 0;
        Cell *finGroupStart {};
        for (int box = 0; box < SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE; box++) {
          Cell *groupStart = &lineIterator[box * SUDOKU_BOX_SIZE];
          auto &groupIterator = (rowMode) ? (HouseIterator &) RowIterator { groupStart }
                                          : (HouseIterator &) ColIterator { groupStart };
          if (candidateCountInGroup( groupIterator, clueIndex ) == 0)
            continue;
          if (killZone_.coversGroup( groupIterator, clueIndex ))
            usualMatch += 1;                               // further fin/sashimi group still possible
          else if (killZone_.absentInChute( &groupIterator[0] )) {
            usualMatch = 0;                                // current line incapable of finned/sashimi N-fish
            break;
          } else                                           // matching (under fin/sashimi conditions)
            finMatch += 1, finGroupStart = groupStart;
        }
        if (usualMatch > 0  &&  finMatch == 1) {           // exactly 1 fin/sashimi per line allowed
          HouseIterator &finGroupIterator = (rowMode) ? (HouseIterator &) RowIterator { finGroupStart }
                                                      : (HouseIterator &) ColIterator { finGroupStart };
          if (finnedFishCandidateCrossOut( finGroupIterator, fishSequence_, killZone_, clueIndex )) {
            printf( "%s with ONE finned/Sashimi group, clue:%d %s: ",
                    (N == 2) ? "x-Wing"    :
                    (N == 3) ? "swordfish" :
                    (N == 4) ? "jellyfish" :
                    (N == 5) ? "starfish"  :
                    (N == 6) ? "sixfish"   :
                    (N == 7) ? "sevenfish" :
                    (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
            for (size_t i = 0; i < N - Step; i++)
              printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
            printf( " (sizes:" );
            for (size_t i = 0; i < N - Step; i++)
              printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
            printf( ")\n  Fin/Sashimi found in line:%d, usualMatch:%d, finGroupStart @[row:%d, col:%d]\n",
                    line, usualMatch, rowNumber( &finGroupIterator[0] ), colNumber( &finGroupIterator[0] ));
            found = true;
          }
        }
      #endif
      }
      return found;
    }  // ---------------------------------------------------------------------------------------

    template<size_t N, size_t Step>
    bool analysisStep( HouseIterator &houseIterator, bool rowMode, int clueIndex, int previousRow ) {
      bool found = false;

      static_assert( Step <= N - 1  &&  Step > 0 );
      //printf( "fishAnalysisStep<%zd, %zd>(mode: %s, clueIndex %d)\n", N, Step, (rowMode) ? "rows" : "cols", clueIndex + 1 );
      for (int tailRow = previousRow + 1; tailRow < elementsof( mergedSum_ ) /* && found == false */; tailRow++) {
        // reduced loop: "tailRow < elementsof( mergedSum_ ) - (Step - 1)" desired for search without fin/sashimi
        // full range:  (without " - (Step - 1)" reasonable for N-fish search with fin/sashimi cells
        if (CellFishSequence<N>::invalidClueCount( mergedSum_[tailRow].possibilities_[clueIndex].count_ ))
          continue;
        fishSequence_[N - Step].CellFishSequenceInit( /* houseIterator, */ rowMode, mergedSum_[tailRow], clueIndex);
        if (killZone_ != fishSequence_[N - Step])          // constellation invalid for N-fish?
          continue;
        killZone_     += fishSequence_[N - Step];
        fishRow_[N - Step] = tailRow;
        //printf( "--> fishAnalysisStep<>(%zd, %zd) recursive call, clue: %d\n", N, Step - 1, clueIndex + 1 );
        if (analysisStep<N, Step - 1>( houseIterator, rowMode, clueIndex, tailRow )) {
          found = true;                // mergedSum_ becomes invalid for current clueIndex after candidate cross out
          break;
        }
      }
    #if 1
      if (Step == 2  &&  found == false)
        if (completeFishSequenceWithTwoFins<N, Step>( rowMode, clueIndex ))
          found = true;
      if (Step == 1  &&  found == false)
        if (completeFishSequenceWithOneFin<N, Step>( rowMode, clueIndex ))
          found = true;
    #endif
      return found;
    }  // ---------------------------------------------------------------------------------------

    template<>
    bool analysisStep<N, 0>( HouseIterator &houseIterator, bool rowMode, int clueIndex, int previousRow ) {
      // This template specialization lets the compiler see, that template's recursion will
      // stop when <Step> reaches a particular value (here: zero).
      // It avoids the C1202 (MS Visual Studio) compiler error: "recursive type or dependency too complex".
      // 
      // (<Step> template argument runs downwards, so it will reach zero.
      //   If <Step> would be running upwards, then a template specialization for a particular <Step>
      //   value greater than start value should also avoid the C1202 error.)

      // N-fish constellation found! Finally try to cross out.
    #if 0
      printf( "%s(clue=%d) %s:",
              (N == 2) ? "x-Wing"    :
              (N == 3) ? "swordfish" :
              (N == 4) ? "jellyfish" :
              (N == 5) ? "starfish"  :
              (N == 6) ? "sixfish"   :
              (N == 7) ? "sevenfish" :
              (N == 8) ? "octopus"   : "??-fish", clueIndex + 1, (rowMode) ? "rows" : "cols" );
      for (size_t i = 0; i < N; i++)
        printf( "%s%d", (i > 0) ? "," : "", fishRow_[i]);
      printf( " (sizes:" );
      for (size_t i = 0; i < N; i++)
        printf( "%s%d", (i > 0) ? "," : "", mergedSum_[fishRow_[i]].possibilities_[clueIndex].count_ );
      printf( ") -- _LAST_ fishAnalysisStep<%zd, 0>()\n", N );
    #endif
      return FishMethods<N>::fishCandidateCrossOut( fishSequence_, killZone_, clueIndex );
    }  // ---------------------------------------------------------------------------------------

  public:
    bool analysis( HouseIterator &&houseIterator, bool rowMode ) {
      // return false, if neither clue possibility can be crossed out
      // N == 2 : x-Wing algorithm level: intermediate. Description:
      //   http://www.kristanix.com/sudokuepic/sudoku-solving-techniques.php
      //   https://www.sudoku.org.pl/skrzydlica.html
      //   YouTube: "dxSudoku #9 The X-Wing"
      // 
      //  ROW mode:                                            COLUMN mode:   (mode depends on houseIterator class)
      //
      //  processing     head.fish_[1] --- head.fish_[0]       processing     head.fish_[1] --- tail.fish_[1]
      //  direction:      |                 |                  direction:      |                 |
      //  top --> bottom  |                 |                  left --> right  |                 |
      //      row         |                 |                     column       |                 |
      //                 tail.fish_[1] --- fishTail.fish_[0]                  head.fish_[0] --- tail.fish_[0]
      //
      //  x-Wing rectangle corner cells: [fishHead.fish_[0], fishTail.fish[1], fishTail.fish_[0], fishTail.fish_[1]]
      //  A x-Wing constellation consists of always 4 cells, which occupy 2 rows and 2 columns.
      // ----------------------------------------------------------------------------------------
      // N == 3: swordfish algorithm level: advanced. Description:
      //   http://www.kristanix.com/sudokuepic/sudoku-solving-techniques.php
      //   https://www.sudoku.org.pl/miecznik.html
      //   YouTube: "dxSudoku #38 Swordfish Puzzle Solving Technique"
      //
      //  ROW mode:
      //
      //  processing       head.fish_[2] ---- head.fish_[1] ---- head.fish_[0]
      //  direction:        |                  |                  |
      //  top --> bottom   body.fish_[2] ---- body.fish_[1] ---- body.fish_[0]
      //      row           |                  |                  |
      //                   tail.fish_[2] ---- tail.fish_[1] ---- tail.fish_[0]
      //
      //  Three swordfish cell-trios: ( head[0,1,2], body[0,1,2], tail[0,1,2] )
      //  However: each triple have to contain at least 2 cells, to keep the swordfish method working.
      //  A swordfish row mode constellation consists of 6,7,8 or 9 cells, which occupy 3 rows and 2..3 columns.
      // 
      //  COLUMN mode:   (mode depends on houseIterator class)
      //
      //  processing       head.fish_[2] ---- body.fish_[2] ---- tail.fish_[2]
      //  direction:        |                  |                  |
      //  left --> right   head.fish_[1] ---- body.fish_[1] ---- tail.fish_[1]
      //      column        |                  |                  |
      //                   head.fish_[0] ---- body.fish_[0] ---- tail.fish_[0]
      //
      //  A swordfish column mode constellation consists of 6,7,8 or 9 cells, which occupy 3 columns and 2..3 rows.
      // ----------------------------------------------------------------------------------------
      // N == 4: jellyfish algorithm level: advanced. Description:
      //   https://www.sudoku.org.pl/meduza.html
      //   YouTube: "dxSudoku #69 Jellyfish Puzzle Solving Technique"
      //
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
      // ----------------------------------------------------------------------------------------
      // N == 5: starfish algorithm level: very advanced. Description:
      //   starfish goes one step further than jellyfish method (and two steps further than swordfish).
      //   An appliable YT video was not found.
      // 
      // N == 6: sixfish algorithm level: very advanced. Description:
      //   Sixfish goes one step further than starfish method.
      //   An appliable YT video was not found.
      // 
      // N == 7: sevenfish algorithm level: very advanced. Description:
      //   Sevenfish goes one step further than sixfish method.
      //   An appliable YT video was not found.
      // 
      // N == 8: octopus algorithm level: very advanced. Description:
      //   Octopus goes one step further than sevenfish method.
      //   An appliable YT video was not found.
      // 
      // In all N-fish methods same principle should be applied on grid rows and grid columns!
      bool found = false;

      for (auto &item : mergedSum_) {
        item.reset();
        for (int line = 0; line < SUDOKU_GRID_SIZE; line++)
          houseIterator[line].mergeCellPossibilities( &item );
        houseIterator.nextHouse();
      }

      for (int clueIndex = 0; clueIndex < elementsof( mergedSum_[0].possibilities_ ); clueIndex++)
        for (int headRow = 0; headRow < elementsof( mergedSum_ ) /* "- (N - 1)" :: full loop for fin search */; headRow++) {
          if (CellFishSequence<N>::invalidClueCount( mergedSum_[headRow].possibilities_[clueIndex].count_ ))
            continue;
          fishSequence_[0].CellFishSequenceInit( /* houseIterator, */ rowMode, mergedSum_[headRow], clueIndex);

          killZone_   = fishSequence_[0];
          /* dxDebug: */
          //    2. ulepszyæ CellFishSequenceInit(), pozbawiaj¹c go argumentu houseIterator, który jest odtwarzany jako lineIterator
          //      z rowMode i z wskaŸnika na pierwsz¹ komórkê ( "mergedSum[] / possibilities" .possibilities_[clueIndex].cell_)
          //      albo przynajmniej wprowadziæ w CellFishSequenceInit() przed wywo³aniem lookBackForCellWithSamePotentialClue()
          //      weryfikacjê (assert) iloœci komórek z clueIndex w current line (przez candidateCountInHouse( lineIt, clueIndex ) )

          fishRow_[0] = headRow;

          if (analysisStep<N, N - 1>( houseIterator, rowMode, clueIndex, headRow )) {
            found = true;
            break;                     // mergedSum_ becomes invalid for current clueIndex after cross out
          }
        }
      return found;
    }  // ---------------------------------------------------------------------------------------

    friend class SudokuGrid;
  };

  using FishMethodDuo     = FishMethods<2>;
  using FishMethodTrio    = FishMethods<3>;
  using FishMethodQuartet = FishMethods<4>;
  using FishMethodQuintet = FishMethods<5>;
  using FishMethodSextet  = FishMethods<6>;
  using FishMethodSeptet  = FishMethods<7>;
  using FishMethodOctet   = FishMethods<8>;

  //================================================================================================
  // todo!
  // 
  // 'std::iterator<std::random_access_iterator_tag, SudokuGrid::Cell, ptrdiff_t, SudokuGrid::Cell *, SudokuGrid::Cell &>':
  // warning STL4015: The std::iterator class template (used as a base class to provide typedefs) is deprecated in C++17.
  // (The <iterator> header is NOT deprecated.)
  // The C++ Standard has never required user-defined iterators to derive from std::iterator. To fix this warning, 
  // stop deriving from std::iterator and start providing publicly accessible typedefs named 
  // iterator_category, value_type, difference_type, pointer, and reference. 
  // Note that value_type is required to be non-const, even for constant iterators.
  // You can define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING or _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS
  // to acknowledge that you have received this warning.
  //
  //template<Cell> class std::iterator_traits;

  enum class HouseId : int { box   = 0, row      = 0x01, col      = 0x02 };
  static bool validHouse( HouseId id ) { return (id >= HouseId::box && id <= HouseId::col ); }
  enum class NodeId  : int { point = 0, rowGroup = 0x01, colGroup = 0x02 };
  static bool validNodeId( NodeId id ) { return (id >= NodeId::point && id <= NodeId::colGroup ); }

  class HouseIterator {
  public:
    typedef Cell                            value_type;
    typedef ptrdiff_t                       difference_type;
    typedef Cell                           *pointer;
    typedef Cell                           &reference;
    typedef std::random_access_iterator_tag iterator_category;
  protected:
    pointer cell_;
    pointer const base_;

    explicit HouseIterator( pointer cell ) : cell_ { cell }, base_ { cell_ } {}

    //bool operator == ( const HouseIterator &rhs ) const { return cell_ == rhs.cell_; }
    //bool operator != ( const HouseIterator &rhs ) const { return cell_ != rhs.cell_; }
    virtual HouseIterator &operator ++ ()    = 0;                              // pre increment
    virtual HouseIterator &operator -- ()    = 0;                              // pre decrement
  //virtual HouseIterator  operator ++ (int) = 0;                              // post increment: impossible in pure virtual
  //virtual HouseIterator  operator -- (int) = 0;                              // post decrement: impossible in pure virtaul
    virtual       reference operator [](difference_type idx) = 0;              // wr access
    virtual const reference operator [](difference_type idx) const = 0;        // rd access
    virtual void nextHouse( void ) = 0;
    virtual void nextGroup( void ) = 0;
    virtual void nextBoxGroup( void ) = 0;
    virtual HouseId type( void ) const = 0;
    void reset( void ) { cell_ = base_; };
  #if 0
    // Arithmetic operators returning own class object don't work with polymorphic classes,
    // a workaround is possible, but expensive:
    // https://stackoverflow.com/questions/37353387/incremental-operator-overload-in-an-abstract-class-c
    // Hence, this kind of operators must be avoided/get around or HouseIterator class cannot be abstract.
    virtual HouseIterator  operator ++ (int)     { HouseIterator tmp = *this; operator++(); return tmp; }
    virtual HouseIterator  operator -- (int)     { HouseIterator tmp = *this; operator--(); return tmp; }
  #endif
          pointer   operator -> ()               { return cell_; }             // wr
    const pointer   operator -> ()         const { return cell_; }             // rd
          reference operator *  ()               { return *cell_; }            // wr
    const reference operator *  ()         const { return *cell_; }            // rd

    HouseIterator &operator += (int rhs) { cell_ += rhs; return *this; }
    HouseIterator &operator -= (int rhs) { cell_ -= rhs; return *this; }

    bool isRowType() const { return type() == HouseId::row; };
    bool isColType() const { return type() == HouseId::col; };
    bool isBoxType() const { return type() == HouseId::box; };
    bool isLineType()const { return type() != HouseId::box; };

    friend class SudokuGrid;
  };

  class RowIterator : public HouseIterator {
    HouseIterator &operator ++ () override { ++cell_; return *this; }                              // pre increment
    HouseIterator &operator -- () override { --cell_; return *this; }                              // pre decrement
    RowIterator   operator ++ (int) { RowIterator t = *this; ++cell_; return t; }                  // post increment
    RowIterator   operator -- (int) { RowIterator t = *this; --cell_; return t; }                  // post decrement
          reference operator [] (difference_type idx)       override { return cell_[idx]; }        // wr access
    const reference operator [] (difference_type idx) const override { return cell_[idx]; }        // rd access
    void nextHouse( void ) override { cell_ += SUDOKU_GRID_SIZE; }             // cannot be used with ++/-- operators!
    void nextGroup( void ) override { cell_ += SUDOKU_BOX_SIZE; }
    void nextBoxGroup( void ) override { cell_ += SUDOKU_GRID_SIZE; }
    HouseId type( void ) const override { return HouseId::row; }

    //static RowIterator iterator( const Cell *cell ) {
    //  RowIterator rowIterator { cell->rowOfCellsStart_ };
    //  return rowIterator;
    //}

    //~RowIterator() {}
    //RowIterator() : HouseIterator( nullptr ) {}
    RowIterator() = delete;
    explicit RowIterator( pointer cell ) : HouseIterator( cell ) {}

    friend class SudokuGrid;
  };
  class ColIterator : public HouseIterator {
    HouseIterator &operator ++ () override { cell_ += SUDOKU_GRID_SIZE; return *this; }            // pre increment
    HouseIterator &operator -- () override { cell_ -= SUDOKU_GRID_SIZE; return *this; }            // pre decrement
    ColIterator   operator ++ (int) { ColIterator tmp = *this; cell_ += SUDOKU_GRID_SIZE; return tmp; }      // post inc
    ColIterator   operator -- (int) { ColIterator tmp = *this; cell_ -= SUDOKU_GRID_SIZE; return tmp; }      // post dec

    reference operator [] (difference_type idx)       override {                                   // wr access
      return cell_[idx * SUDOKU_GRID_SIZE];
    }
    const
    reference operator [] (difference_type idx) const override {                                   // rd access
      return cell_[idx * SUDOKU_GRID_SIZE];
    }
    void nextHouse( void ) override { cell_ += 1; }                            // cannot be used with ++/-- operators!
    void nextGroup( void ) override { cell_ += SUDOKU_BOX_SIZE * SUDOKU_GRID_SIZE; }
    void nextBoxGroup( void ) override { cell_ += 1; }
    HouseId type( void ) const override { return HouseId::col; }

    //~ColIterator() {}
    //ColIterator() : HouseIterator( nullptr ) {}
    ColIterator() = delete;
    explicit ColIterator( pointer cell ) : HouseIterator( cell ) {}

    friend class SudokuGrid;
  };
  class BoxIterator : public HouseIterator {
    HouseIterator &operator ++ () override {                                                        // pre increment
      ++cell_;
      if ((cell_ - base_) % SUDOKU_BOX_SIZE == 0)          // next row of the current box?
        cell_ += SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE;
      return *this;
    }
    HouseIterator &operator -- () override {                                                        // pre decrement
      cell_ -= ((cell_ - base_) % SUDOKU_BOX_SIZE == 0) ? SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE + 1 : 1;
      return *this; 
    }
    BoxIterator operator ++ (int) { BoxIterator tmp = *this; operator++(); return tmp; }           // post increment
    BoxIterator operator -- (int) { BoxIterator tmp = *this; operator--(); return tmp; }           // post decrement
    reference operator [] (difference_type idx) override {                                         // wr access
    #if 1
      static_assert(SUDOKU_GRID_SIZE == 9);
      static constexpr difference_type offset[SUDOKU_GRID_SIZE] = { 
        0,
        1 + (1 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        2 + (2 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        3 + (3 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        4 + (4 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        5 + (5 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        6 + (6 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        7 + (7 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
        8 + (8 / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE),
      };
      return cell_[offset[idx]];
    #else
      return cell_[idx + (idx / SUDOKU_BOX_SIZE) * (SUDOKU_GRID_SIZE - SUDOKU_BOX_SIZE)];
    #endif
    }
    const reference operator [] (difference_type idx) const override {                             // rd access
      return operator[]( idx );
    }
    void nextHouse( void ) override {                      // cannot be used with ++/-- operators!
      cell_ += SUDOKU_BOX_SIZE;
      if ((cell_ - base_) % SUDOKU_GRID_SIZE == 0)         // next row of boxes?
        cell_ += SUDOKU_GRID_SIZE * (SUDOKU_BOX_SIZE - 1);
    }
    void nextGroup( void ) override { assert( 0 ); }       // illegal, lack of information about group node orientation
    void nextBoxGroup( void ) override { assert( 0 ); }    // illegal, lack of information about group node orientation
    //void nextRowGroup( void ) { cell_ += SUDOKU_GRID_SIZE; }
    //void nextColGroup( void ) { cell_ += 1; }
    //void nextGroup( bool rowMode ) { (rowMode) ? nextRowGroup() : nextColGroup(); }
    HouseId type( void ) const override { return HouseId::box; }

    void nextChuteBox( bool rowMode ) {                    // traverse boxes across/horizontally or down/vertically chutes
      cell_ += (rowMode) ? SUDOKU_BOX_SIZE : SUDOKU_BOX_SIZE * SUDOKU_GRID_SIZE;
    }
    bool operator == ( const Cell *const rhs) const {
      return (this->cell_ == rhs);
    }

    //~BoxIterator() {}
    //BoxIterator() : HouseIterator( nullptr ) {}
    BoxIterator() = delete;
    explicit BoxIterator( pointer cell ) : HouseIterator( cell ) {}

    friend class SudokuGrid;
  };
#if 0
  static HouseIterator &lineIterator( bool rowMode, const Cell *cell ) {
    --> problem: zwracana jest referencja na obiekt ¿yj¹cy tylko w czasie wywo³ania lineIterator()
    //          jak napisaæ kod, aby zwracaæ wybrany iterator i zainicjalizowaæ virtualny HouseIterator w kodzie wywo³uj¹cym??
    return (rowMode) ? (HouseIterator &) RowIterator { const_cast<Cell *>( cell ) }
                     : (HouseIterator &) ColIterator { const_cast<Cell *>( cell ) };
#elif 0
  --> czym poni¿szy wariant siê ró¿ni?
  static HouseIterator &&lineIterator( bool rowMode, const Cell *cell ) {
    // problem: zwracana jest referencja na rvalue 
    return (rowMode) ? (HouseIterator &&) RowIterator { const_cast<Cell *>( cell ) }
                     : (HouseIterator &&) ColIterator { const_cast<Cell *>( cell ) };
  }
#endif

  const char *const riddleName_;
  const uchar (*const clues_)[SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE];

  std::array<Cell, SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE> *grid_;

  // todo! usun¹æ niepotrzebne warianty chainProtection(), pozosta³e PRZENIEŒÆ do SudokuSolver.cpp
#if 1 // required by disjoinChainAnalysis()
  static void chainProtection( bool state, Cell &head, Cell &body, Cell &tail ) {
    head.candidateProtection( state );
    body.candidateProtection( state );
    tail.candidateProtection( state );
  }
#endif

  static void chainProtection( bool state, Cell **cells, size_t size ) {
    for (; size-- > 0; (*cells++)->candidateProtection( state ));
  }
  static void candidateHouseProtection( bool state, HouseIterator &houseIterator );
  static void candidateGroupProtection( bool state, HouseIterator &groupIterator );
  static bool candidateGroupProtected  ( HouseIterator &groupIterator );
  static bool candidateGroupUnprotected( HouseIterator &groupIterator );
  static bool candidateGroupUnprotected( HouseIterator &groupIterator, int clueIndex );
  static bool candidateCrossOutInProtectedHouse( HouseIterator &houseIterator, int clueIndex );
  static bool candidateCrossOutInProtectedHouse( HouseIterator &&houseIterator, int clueIndex ) {
    return (candidateCrossOutInProtectedHouse( houseIterator, clueIndex ));
  }
  static bool candidateCrossOutInProtectedHouse( HouseIterator &houseIterator, const Cell::CandidateTraits &candidateTraits );
  static bool candidateCrossOutInProtectedHouse( HouseIterator &&houseIterator, const Cell::CandidateTraits &candidateTraits ) {
    return (candidateCrossOutInProtectedHouse( houseIterator, candidateTraits ));
  }
  static bool candidateCrossOutInProtectedHouseIfInAllNeighboursSight( HouseIterator &houseIterator, int clueIndex,
                                                                       const Cell *const *neighbours, size_t neighboursQuantity );
  static bool candidateCrossOutInProtectedHouseIfInAllNeighboursSight( HouseIterator &&houseIterator, int clueIndex,
                                                                       const Cell *const *neighbours, size_t neighboursQuantity ) {
    return (candidateCrossOutInProtectedHouseIfInAllNeighboursSight( houseIterator, clueIndex, neighbours, neighboursQuantity ));
  }
  static void candidateCrossOutInRowColBox( Cell *cell );
  static void establishClueAndCrossItOutInRowColBox( Cell *cell, int clue );

  static int  candidateCountInHouse( HouseIterator &houseIterator, int clueIndex );
  static int  candidateCountInHouse( HouseIterator &&houseIterator, int clueIndex ) {
    return candidateCountInHouse( houseIterator, clueIndex );
  }
  static int  candidateCountInGroup( HouseIterator &groupIterator, int clueIndex );
  static int  candidateCountInGroup( HouseIterator &&groupIterator, int clueIndex ) {
    return candidateCountInGroup( groupIterator, clueIndex );
  }
  static std::tuple<int, int> candidateCountInGroup( HouseIterator &groupIterator1,
                                                     HouseIterator &groupIterator2, int clueIndex ) {
    return { candidateCountInGroup( groupIterator1, clueIndex ),
             candidateCountInGroup( groupIterator2, clueIndex ) };   // RVO - cheers!
  }

  void initGrid( void );

  bool singleSolution( void );

  static bool singleCellAnalysis( HouseIterator &&houseIterator );
  bool singleCell( void );

  void singleBoxAnalysisInit( HouseIterator &houseIterator,
                              Possibilities  (*mergeSumData)[SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE] );
  bool singleBoxAnalysis    ( HouseIterator &houseIterator,
                              Possibilities  (*mergeSumData)[SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE],
                              Possibilities *(*mergeSum)    [SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE] );
  bool singleBox( void );

  bool disjoinSubsetAnalysis( HouseIterator &houseIterator );
  bool disjoinSubset( void );

  bool disjoinChainAnalysis( HouseIterator &houseIterator );
  bool disjoinChain( void );

  //bool nakedDuoAnalysis( HouseIterator &&houseIterator );
  bool nakedDuo( void );
  //bool nakedTrioAnalysis( HouseIterator &&houseIterator );
  bool nakedTrio( void );
  bool nakedQuartet( void );
  bool nakedQuintet( void );
  bool nakedSextet( void );
  bool nakedSeptet( void );
  bool nakedOctet( void );

  static Cell *lookBackForCellWithSamePotentialClue( HouseIterator &lineIterator, int clueIndex );
  static Cell *lookBackForCellWithSamePotentialClue( HouseIterator &houseIterator, Cell *seedCell, int clueIndex );

  bool x_Wing( void );
  bool swordfish( void );
#if 0
  bool jellyfishCrossOut( CellFishQuartet (&fishQuartet)[4], int clueIndex);   // without recursive templated method
  bool jellyfishAnalysis( HouseIterator &&houseIterator, bool rowMode );       // without recursive templated method
#endif
  bool jellyfish( void );
#if 0
  bool starfishCrossOut( CellFishQuintet (&fishQuintet)[5], int clueIndex);    // without recursive templated method
  bool starfishAnalysis( HouseIterator &&houseIterator, bool rowMode );        // without recursive templated method
#endif
  bool starfish( void );     // aka: squirmbag (source: YT: Sudoku Swami tutorial #10, mark 3:30)
  bool sixfish( void );      // aka: whale     (source: YT: Sudoku Swami tutorial #10, mark 3:30)
  bool sevenfish( void );    // aka: leviathan (source: YT: Sudoku Swami tutorial #10, mark 3:30)
  bool octopus( void );

  static bool w_WingCrossOut( const Cell *firstEP, const Cell *otherEP, int killCandidateIndex );
  static bool w_WingTraversing( RowIterator &&rowIterator, Cell *(*quartet)[4], int (*candidateXY)[2] );
  static bool w_WingTraversing( RowIterator &&rowIterator, Cell *(*quartet)[4] );
  static bool w_WingAnalysis( RowIterator &&rowIterator );
  bool w_Wing( void );

  static bool xy_WingCrossOut( const Cell *pivot, const Cell *firstPincer, const Cell *otherPincer, int pincerZ );
  static bool xy_WingAnalysis( HouseIterator &otherPincerIterator, const Cell *pivot, const int (*pivotXY)[2], const Cell *firstPincer );
  static bool xy_WingAnalysis( RowIterator &&rowIterator );
  bool xy_Wing( void );

  static bool xyz_WingAnalysis( HouseIterator &&otherPincerIterator, Cell *(*trio)[3], Cell::CandidateTraits &otherPincerForm );
  static bool xyz_WingAnalysis( RowIterator &&rowIterator );
  bool xyz_Wing( void );

  static bool skyscraperInvalidConstellation( const CellFishDuo &pair );
  static bool skyscraperInvalidConstellation(const CellFishDuo &head, const CellFishDuo &tail);
  static bool skyscraperCrossOut( Cell *head, Cell *tail, bool rowMode, int clueIndex );
  static bool skyscraperCrossOut( CellFishDuo *head, CellFishDuo *tail, int clueIndex );
  static bool skyscraperAnalysis( HouseIterator &&houseIterator, bool rowMode );
  bool skyscraper( void );             // == variant of X-Chain

  static bool twoStringKiteCrossOut( const Cell *firstEP, const Cell *otherEP, int killCandidateIndex );
  static bool twoStringKiteInvalidConstellation( const CellFishDuo &xWingHalf );
  static bool twoStringKiteInvalidConstellation( const CellFishDuo &xWingHead, const CellFishDuo &xWingTail );
  static bool twoStringKiteAnalysis( ColIterator &&colIterator, const CellFishDuo &xWingHead, int clueIndex );
  static bool twoStringKiteAnalysis( RowIterator &&rowIterator );
  bool twoStringKite( void );          // == variant of X-Chain

  static bool emptyRectangleInvalidConstellation( const CellFishDuo &pair );
  static Cell *emptyRectangleKillZone( BoxIterator *boxOfCells, const CellFishDuo &head, int clueIndex );
  static bool emptyRectangleCrossOut( Cell *(*killZone)[2 * (SUDOKU_GRID_SIZE / SUDOKU_BOX_SIZE - 1)], int clueIndex );
  static bool emptyRectangleAnalysis( HouseIterator &&houseIterator, bool rowMode );
  bool emptyRectangle( void );         // == variant of X-Chain

  static bool finnedX_WingAnalysis( const CellFishDuo &xWingHead, const Cell *third, bool rowMode, int clueIndex );
  static bool finnedX_WingAnalysis( HouseIterator &&houseIterator, bool rowMode );
  bool finnedX_Wing( void );           // == variant of X-Chain
  //bool sashimiX_Wing( void );        // covered by Finned X-Wing implementation, see: dxsudoku #76 + Practice Exercises #11
  bool finnedSwordfish( void );
  bool finnedFish( void );

  static bool uniqueRectangleAnalysisType1( HouseIterator &houseIterator,  CellFishQuartet &rectangle );
  static bool uniqueRectangleAnalysisType1( HouseIterator &&houseIterator, CellFishQuartet &rectangle ) {
    return (uniqueRectangleAnalysisType1( houseIterator, rectangle ));
  }  // -------------------------------------------------------------------------------------------
  static bool uniqueRectangleAnalysisType2( HouseIterator &houseIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleAnalysisType3( HouseIterator &houseIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleAnalysisType4( HouseIterator &houseIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleAnalysis( HouseIterator &&houseIterator, bool rowMode );
  static bool uniqueRectangleType5CrossOut( const Cell *rectangleCorner, CellFishQuartet &rectangle, int clueIndex );
  static bool uniqueRectangleScrutinyType5( RowIterator &rowIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleScrutinyType6( RowIterator &rowIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleScrutinyType7( RowIterator &rowIterator, CellFishQuartet &rectangle );
  static bool uniqueRectangleScrutiny( RowIterator &&rowIterator, bool rowMode );
  bool uniqueRectangle( void );

#if 0
  bool hiddenRectangle( void );        // == Unique Rectangle Type 7 (YT: "dxSudoku #56 Hidden Rectangle")

  bool turbotFish( void );             // == variant of X-Chain with exactly 4 cells and 3 links
#endif

  bool x_Chain( void );                // dxSudoku: #41, #73, #100, #94, #93

  class X_Chain {
    class Link {
      const Cell *cell_;
      const Link *const prev_;
      static bool chainKillZone( X_Chain *chain, const Cell *endpoint, const Cell *strongLinkStart, HouseIterator &&houseIterator );
      static bool chainKillZone( X_Chain *chain, const Cell *endpoint, const Cell *strongLinkStart );
             bool chainKillZone( X_Chain *chain, const Cell *endpoint );
      static bool chainClosedKillZone( X_Chain *chain, const Cell *endpoint, const Cell *strongLinkStart );
             bool chainClosedKillZone( X_Chain *chain /* , const Cell *endpoint */);

      bool isStrongLinkStartInChain( const Cell *cell );
      bool isStrongLinkEndInChain  ( const Cell *cell ) { return isStrongLinkStartInChain( cell ) == false; }
      void printChain( const char *title, const Cell *cell, const Link *link, const char *msg );
      bool addStrongLink( X_Chain *chain, HouseIterator &&houseIterator );
      bool addStrongLink( X_Chain *chain );
      bool addWeakLink  ( X_Chain *chain, HouseIterator &&houseIterator );
      bool addWeakLink  ( X_Chain *chain );

      Link() = delete;
      Link( const Cell *cell, const Link *previous ) : cell_ { cell }, prev_ { previous } {}

    //friend bool             X_Chain::seedChain( HouseIterator &&houseIterator, int clueIndex );
    //friend bool SudokuGrid::X_Chain::seedChain( HouseIterator &&houseIterator, int clueIndex );  full name not required
      friend class X_Chain;
    };

    int   clueIndex_;
    Link  seedLink[2];                 // seedLink[1].cell_ == chain's root cell
    //size_t      killZoneCnt;
    //Cell       *killZone[20];

    void  setClueIndex( int clueIndex ) { clueIndex_ = clueIndex; }
    int   getClueIndex( void ) const { return clueIndex_; }
    const Cell *getRoot( void ) const { return seedLink[1].cell_; }
    bool  seedStrongLink( HouseIterator &houseIterator /*, Link(*seedLink)[2] */ );

    X_Chain() = delete;
    X_Chain( int clueIndex ) :
      clueIndex_ { clueIndex },
      seedLink {{ nullptr, &seedLink[1] }, { nullptr, nullptr }} {}

  public:
    static bool seedChain( HouseIterator &&houseIterator, int clueIndex );
  //friend bool SudokuGrid::x_Chain( void );
  };

  bool x_NodeChain( void );            // dxSudoku: #55: 10:30, #75: 5:35, #77: 1:30 (Sashimi chain may be easily extended!)

  class X_NodeChain {
    class Link {
      struct Node {
        NodeId     type_;
        const Cell *cell_;

        void candidateProtection( bool state );
        bool sharingSameRow( const Node *node ) const;
        bool sharingSameCol( const Node *node ) const;
        bool sharingSameBox( const Node *node ) const;
        bool sharingSameBox_sameEntityAllowed( const Node *node ) const;
        bool sharingSameHouse( const Node *node ) const;
        bool sharingSameHouse( const Cell *cell ) const;
        bool occupyingDifferentRow( const Node *node ) const { return sharingSameRow( node ) == false; }
        bool occupyingDifferentCol( const Node *node ) const { return sharingSameCol( node ) == false; }
        bool occupyingDifferentBox( const Node *node ) const { return sharingSameBox( node ) == false; }
        bool occupyingDifferentBox_sameEntityAllowed( const Node *node ) const {
          return sharingSameBox_sameEntityAllowed( node ) == false;
        }
        bool isPointNode() const { return type_ == NodeId::point; }
        bool isGroupNode() const { return isPointNode() == false; }
        bool isRowGroupNode() const { return type_ == NodeId::rowGroup; }
        bool isColGroupNode() const { return type_ == NodeId::colGroup; }
        bool isPerpendicular( const NodeId &nodeId ) const {
          assert( type_ != NodeId::point && nodeId != NodeId::point );
          static_assert(0x01 == static_cast<int>( NodeId::rowGroup )  &&
                        0x02 == static_cast<int>( NodeId::colGroup ));
          return (static_cast<int>( type_ ) | static_cast<int>( nodeId )) ==
                 (static_cast<int>( NodeId::rowGroup ) | static_cast<int>( NodeId::colGroup ));
        }
        bool isPerpendicular( const Node &node ) const { return isPerpendicular( node.type_ ); }
        bool isPerpendicular( const HouseId houseId ) const {
          return isPerpendicular( static_cast<NodeId>( houseId ) );
        }
        bool isParallel( const NodeId  nodeId  ) const { return isPerpendicular( nodeId  ) == false; }
        bool isParallel( const Node   &node    ) const { return isPerpendicular( node    ) == false; }
        bool isParallel( const HouseId houseId ) const { return isPerpendicular( houseId ) == false; }
        const Cell *intersectionOfNodes( HouseIterator &lineIterator ) const {
          assert( isGroupNode() &&  isPerpendicular( lineIterator.type() ) );
          assert( cell_->boxOfCellsStart_ == lineIterator->boxOfCellsStart_ );           // avoid sharingSameBox() call

          return (isRowGroupNode()) ? intersectionOfCells( cell_, &lineIterator[0] )
                                    : intersectionOfCells( &lineIterator[0], cell_ );
        }
        const Cell *intersectionOfNodes( const Node &node ) const {
          assert( isGroupNode() && node.isGroupNode() && isPerpendicular( node ) );
          assert( cell_->boxOfCellsStart_ == node.cell_->boxOfCellsStart_ );

          return (isRowGroupNode()) ? intersectionOfCells( cell_, node.cell_ )
                                    : intersectionOfCells( node.cell_, cell_ );
        }
        bool isPartOfNode( const Cell *cell ) const;
        void print( const char *prefix_format, ... ) const;
        void printTitled( const char *prefix, const char *suffix_format, ... ) const;
      } node;
      const Link *const prev_;

      static int candidateCountInNode( const Node &node, int clueIndex ) {
        if (node.isPointNode())
          return 1;
        else {
          auto &lineIterator = (node.isRowGroupNode()) ? (HouseIterator &) RowIterator { const_cast<Cell *>( node.cell_ ) }
                                                       : (HouseIterator &) ColIterator { const_cast<Cell *>( node.cell_ ) };
          return candidateCountInGroup( lineIterator, clueIndex );
        }
      }
      size_t length( void ) const;
      void printChainEngine( X_NodeChain *chain,             const char *format, va_list &args ) const;
      void printChain( X_NodeChain *chain,                   const char *format, ... ) const;
      void printChain( X_NodeChain *chain, const Node &node, const char *format, ... ) const;
      static const Cell *pointFromGroup( const Cell *groupStart, bool rowType, int clueIndex, const Cell *cellToIgnore = nullptr );
      bool diveDeeper( X_NodeChain *chain, const Node &foundNode ) const;
      bool chainClosedKillZone( X_NodeChain *chain, const Node &loopStartChainedNode ) const;
      bool chainClosedKillZone( X_NodeChain *chain, const Node &loopStartChainedNode, const Node &endNode ) const;
      static bool killZone( X_NodeChain *chain, const Node &nodeA, const Node &nodeB, HouseIterator &&houseIt);
      static bool killZone( X_NodeChain *chain, const Node &nodeA, const Node &nodeB );
      bool chainKillZone( X_NodeChain *chain, const Node &endNode ) const;
      const Node *getChainedNode( const Cell *cell ) const;
      const Link *chainEffectivelyClosed( const Node &strongLinkEnd ) const;
      bool isStrongLinkStartInChain( const Node &node ) const;
      bool usualEliminationAndDiveDeeper( X_NodeChain *chain, Node &foundNode ) const;
      bool specialElimination( X_NodeChain *chain, Node &foundNode ) const;
      bool eliminate( X_NodeChain *chain, Node &foundNode ) const;
      bool addStrongLinkWithinLine( X_NodeChain *chain, HouseIterator &lineIterator ) const;
      bool addStrongLinkWithinBox( X_NodeChain *chain, HouseIterator &&boxGroupIterator ) const;
      bool addStrongLink( X_NodeChain *chain, HouseIterator &houseIterator ) const;
      bool addStrongLink( X_NodeChain *chain, HouseIterator &&houseIterator ) const { return addStrongLink( chain, houseIterator ); }
      bool addStrongLink( X_NodeChain *chain ) const;
      bool addGroupedNodeWithinLineWeakLink( X_NodeChain *chain, HouseIterator &lineIterator ) const;
      bool addGroupedNodeWithinBoxWeakLink( X_NodeChain *chain, HouseIterator &&boxGroupIterator ) const;
      bool addWeakLink  ( X_NodeChain *chain, HouseIterator &houseIterator ) const;
      bool addWeakLink  ( X_NodeChain *chain, HouseIterator &&houseIterator ) const { return addWeakLink( chain, houseIterator ); };
      bool addWeakLink  ( X_NodeChain *chain ) const;
      void candidateProtection( bool state ) { node.candidateProtection( state ); };
      void arrangeSeedLink( HouseIterator &tailGroupIt, int tailCandidateCount,
                            HouseIterator &headGroupIt, int headCandidateCount,
                            int clueIndex, const Cell *intersectionCell );
      Link() = delete;
      Link( const Cell *cell, NodeId type, const Link *previous ) : node { type, cell }, prev_ { previous} {}

      //static Node getNode( HouseIterator &groupIterator, int candidateCount, int clueIndex );

      //friend bool             X_NodeChain::seedChain( HouseIterator &&houseIterator, int clueIndex );
      //friend bool SudokuGrid::X_NodeChain::seedChain( HouseIterator &&houseIterator, int clueIndex );  full name not required
      friend class X_NodeChain;
    };
    int  clueIndex_;
    Link seedLink[2];        // seedLink[1].node == chain's root

    void  setClueIndex( int clueIndex ) { clueIndex_ = clueIndex; }
    int   getClueIndex( void ) const { return clueIndex_; }

    static Link::Node getNode( HouseIterator &groupIterator, int candidateCount, int clueIndex, const Cell *intersectionCell );
    static Link::Node getNode( HouseIterator &groupIterator, int candidateCount, int clueIndex );

    const Link *getRoot( void ) const { return &seedLink[1]; }

    static void printSeedLink( const Link (&seedLink)[2], const char *format, ... );

    bool  cultureChain( Link (&seedLink)[2] );
    bool  seedGroupStrongLinkPerpendicular( BoxIterator &boxIterator, int candidateCount );
    bool  seedGroupStrongLinkParallel( HouseIterator &&groupIterator );
    bool  seedGroupStrongLinkParallel( BoxIterator &boxIterator );
    bool  seedGroupStrongLink( BoxIterator &boxIterator );
    bool  seedGroupStrongLink( HouseIterator &lineIterator );
    bool  seedStrongLink( HouseIterator &houseIterator );

    X_NodeChain() = delete;
    X_NodeChain( int clueIndex ) :
      clueIndex_ { clueIndex },
      seedLink {{ nullptr, NodeId::point, &seedLink[1] }, { nullptr, NodeId::point, nullptr }} {}

  public:
    static bool seedChain( HouseIterator &&houseIterator, int clueIndex );
  };

  //bool AICs( void );                 // dxSudoku #55: 12:00, #96: 17:20, (AIC = Alternate Inference Chain)
  //bool niceLoop_AIC( void );         // dxSudoku #97, #98, #99, #94
  //bool groupNiceLoop_AIC( void );
  //bool xy_Chain( Void );             // dxSudoku #66 + #67
  //bool remotePair( void );
  //bool sueDeCoq( void );             // possible are 4,5,6,7,8-cell types - (dxSudoku #101 = 4-cell type)
  //bool bowmansBingo( void );         // aka?: "Bivalue Elimination Part I" - dxSudoku #83, (Part II: #84)
  //bool almostLockedSets( void );
  //bool forcingChain( void );
  //bool squirmbag( void );            // aka: starfish
  //bool nishio( void );               // aka: ??

  static bool validClue( Cell &cell, uchar clue );
  bool bruteForce( size_t startIndex );
  bool bruteForce( void );

  static bool sortedBruteForce( std::multimap<int, Cell *> &grid,
                                std::multimap<int, Cell *>::iterator startIterator );
  bool sortedBruteForce( void );

  void printSudokuGrid( const char *format, ... ) const;
  bool sudokuFinished( void );
  bool sudokuSolvingCheck( void );

  static int rowNumber( const Cell *cell ) {
    return (int) (cell - cell->colOfCellsStart_) / SUDOKU_GRID_SIZE;
  }  // -------------------------------------------------------------------------------------------
  static int colNumber( const Cell *cell ) {
    return (int)(cell - cell->rowOfCellsStart_);
  }  // -------------------------------------------------------------------------------------------
  static int boxNumber( const Cell *cell ) {
  #if 1
    Cell *boxOrigin = cell->boxOfCellsStart_;
    int  boxRow     = (int) (boxOrigin - boxOrigin->colOfCellsStart_) / SUDOKU_GRID_SIZE;
    int  boxCol     = (int) (boxOrigin - boxOrigin->rowOfCellsStart_) / SUDOKU_BOX_SIZE;
    return boxRow + boxCol;
  #else
    return (int) (cell->boxOfCellsStart_ - &(*grid_)[0]) / SUDOKU_GRID_SIZE +
           (int) (cell - &(*grid_)[0]) % (SUDOKU_GRID_SIZE * SUDOKU_BOX_SIZE) / SUDOKU_BOX_SIZE;
  #endif
  }  // -------------------------------------------------------------------------------------------

  static Cell *intersectionOfCells( const Cell *rowSource, const Cell *colSource ) {
    return rowSource->rowOfCellsStart_ + (colSource - colSource->rowOfCellsStart_);
  }  // -------------------------------------------------------------------------------------------

  static Cell *intersectionOfCells( const Cell *first, const Cell *second, bool rowMode ) {
    if (rowMode)
      return (first->rowOfCellsStart_  + (second - second->rowOfCellsStart_));
    else
      return (second->rowOfCellsStart_ + (first  -  first->rowOfCellsStart_));
  }  // -------------------------------------------------------------------------------------------

  static std::tuple<const Cell *, bool> intersection( const Cell *rowSource, const Cell *colSource, int clueIndex ) {
    Cell *intersectionCell = intersectionOfCells( rowSource, colSource );
    return { intersectionCell, intersectionCell->candidatePossible( clueIndex ) };
  }  // -------------------------------------------------------------------------------------------

  static Cell *groupNodeStartCell( const Cell *cell, bool rowMode ) {
    return intersectionOfCells( cell, cell->boxOfCellsStart_, rowMode );
  }  // -------------------------------------------------------------------------------------------
#if 0
  static Cell (*groupOfCells( Cell *group ))[SUDOKU_BOX_SIZE] {      // static auto groupOfCells( Cell *group ) {....}
    return reinterpret_cast<Cell (*)[SUDOKU_BOX_SIZE]> (group);
  }  // -------------------------------------------------------------------------------------------
#endif
  static Cell (*groupOfCells( const Cell *cell, bool rowMode ))[SUDOKU_BOX_SIZE] {
    // return pointer to group of 3 cells address of first cell in box row/column (depending on rowMode)
    Cell *group = (rowMode) ? cell->rowOfCellsStart_ + (cell->boxOfCellsStart_ - cell->boxOfCellsStart_->rowOfCellsStart_)
                            : cell->boxOfCellsStart_ + (cell->colOfCellsStart_ - cell->boxOfCellsStart_->colOfCellsStart_);
    return reinterpret_cast<Cell (*)[SUDOKU_BOX_SIZE]> (group);
  }  // -------------------------------------------------------------------------------------------

public:
  ~SudokuGrid();
  SudokuGrid() = delete;                                             // no default constructor
  SudokuGrid( const char *name, const uchar (*input)[SUDOKU_GRID_SIZE * SUDOKU_GRID_SIZE] );

  SudokuGrid( const SudokuGrid &rhs ) = delete;                      // no copy constructor
  SudokuGrid(       SudokuGrid &&rhs ) = delete;                     // no move constructor
  SudokuGrid &operator = ( const SudokuGrid &rhs ) = delete;         // no copy assignment
  SudokuGrid &operator = (       SudokuGrid &&rhs ) = delete;        // no move assignment

  int sudokuSolve( int riddleId = -1 );
};

#endif  // !_SUDOKU_SOLVER_H
