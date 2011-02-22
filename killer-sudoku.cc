/*
  "Killer sudoku" solver. For details of the game, see e.g.
  http://www.timesonline.co.uk/article/0,,7-1757275,00.html
  
  For instructions on using the program, see the attached README file and
  example inputs.

  Some definitions are in order: GRID: the 9x9 field which has to be
  filled with. It comprises 9*9=81 CELLS.  The grid is also subdivided
  in 3x3 SUBGRIDS. It should be obvious what are ROWS and COLUMNS. A
  REGION is a connected set of cells which is not neccessarily
  rectangular.  The sum of the numbers in the region must equal the
  sum specified by the puzzle (TARGET SUM).

  Cells with a value that has not yet been determined are set to 0.
 
  =============================================================================
  Copyright (c) 2005 Zeljko Vrba <zvrba@globalnet.hr>

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <bitset>
#include <stdexcept>

using namespace std;

/*
  The region shape is irrelevant, since addition is associative. Another
  observation: no matter how many 0 cells are there in a region, there is
  always a lower bound at which the program starts its guessing. This lower
  bound can be larger than 1. An example: given a region of 2 empty cells
  with target sum of 11, the lower bound is 2. For if the program tries to
  put 1 in a cell, the sum of 11 can't be achieved since the max. number
  that can be put in the remaining cell is 9, giving a maximum sum of 10.

  Since the shape is irrelevant, the data in this structure is enough to
  keep track of the state during the search for solution.
*/
struct region {
	unsigned	size_;			// number of cells in the region
	unsigned	target_sum_;	// target sum of cells in the region
	unsigned	current_sum_;	// current sum of the cells in the region
	unsigned	nzero_;			// number of zero cells in the region

	region(unsigned target_sum, unsigned size) : 
		size_(size), target_sum_(target_sum), current_sum_(0), nzero_(size)
	{ }

	void add(int n);
	unsigned lowest_feasible_element() const;
};

/*
  Add or subtract n from the current sum in the region. If n is negative,
  then increase the number of zero cells, otherwise decrease it. Also
  performs some sanity checks. n is allowed to be 0, in which case nothing
  is done.
*/
void region::add(int n)
{
	if(n) {
		current_sum_ += n;
		nzero_ += n > 0 ? -1 : 1;
	}

	/* consistency check */
	if((current_sum_ > target_sum_) || (nzero_ > size_))
		throw logic_error("region constraints violated");
}

/*
  As explained above, each region has a lowest element that is feasible to be
  put in any cell. This lower bound, l, is feasible if

  current_sum + l + (nzero-1)*9 >= target_sum

  The multyplying factor is 9 since the numbers in the region MAY repeat. If
  l is negative, just return 1, since all cells are in range [1,9].
*/
unsigned region::lowest_feasible_element() const
{
	int l = target_sum_ - current_sum_ - (nzero_-1)*9;
	return l <= 0 ? 1 : l;
}

/*
  The grid has dimensions 9x9, a total of 81 cells. Each cell is either 0
  (unfilled, to be guessed by the program) or a number 1-9. For each cell
  the grid also maintains a pointer to the region containing the cell,
  allowing fast updates of the region sum.

  Row, column and subgrid validity tests (i.e. that each contains each
  number 1-9) are done with the help of bitsets. For each row, column
  and subgrid, a bitset of 10 elements is maintained (0th element is
  unused). Whenever a number is present in the row/col/subg, the
  corresponding bit in the bitvector is set.

  Cells and subgrids are addressed with a single index, thus
  neccessitating conversion to row and columns in few cases.
*/
struct grid {
	vector<unsigned> cells_;						// all cells
	bitset<10>		rows_[9], columns_[9], sub_[9];	// quick validity test
	vector<region*>	regions_;						// cell -> region map

	grid() : cells_(81, 0), regions_(81, (region*)0)
	{ }

	bool put(unsigned, unsigned);
	void add_region(unsigned, const vector<unsigned>&);
};

/*
  Put number n into cell indexed by idx. If no constraints are violated, the
  method puts the number in the cell (updating its corresponding region) and
  returns true, otherwise it returns false. Putting 0 clears the cell.
*/
bool grid::put(unsigned idx, unsigned n)
{
	/* calculate the index of the 3x3 subgrid corresponding to idx */
	unsigned row = idx / 9;				// row and column in the grid
	unsigned col = idx % 9;
	unsigned subrow = idx / 27;			// row and column in the 3x3 subgrid
	unsigned subcol = (idx % 9) / 3;	// corresponding to idx
	unsigned subidx = 3*subrow+subcol;	// index in the subgrid
	unsigned oldn = cells_[idx]; 		// existing number in the cell
	region *r = regions_[idx];			// region corresponding to the cell

	/* consistency check */
	if((idx > 80) || (n > 9))
		throw invalid_argument("invalid arguments to grid::put");

	/* the specified number already exists in row, col or 3x3 subgrid */
	if(rows_[col][n] || columns_[row][n] || sub_[subidx][n])
		return false;

	/* adding the new number would exceed the target sum */
	if(r->current_sum_ - oldn + n > r->target_sum_)
		return false;

	/* remove existing number from the grid */
	cells_[idx] = 0; r->add(-oldn);
	rows_[col][oldn] = columns_[row][oldn] = sub_[subidx][oldn]  = false;

	/* put new number in the grid */
	if(n > 0) {
		cells_[idx] = n; r->add(n);
		rows_[col][n] = columns_[row][n] = sub_[subidx][n]  = true;
	}

	return true;
}

/*
  Add a new region to the grid. The region is described by the target sum
  and a list of indices in the grid contained in the region.
*/
void grid::add_region(unsigned target_sum, const vector<unsigned> &cidx)
{
	/* create region and update back-references */
	region *r = new region(target_sum, cidx.size());
	for(unsigned i = 0; i < cidx.size(); i++)
		regions_[cidx[i]] = r;
}

/*
  Load data, performing sanity checking along the way. The data format is
  described in the README file. The method returns true on success, and false
  on any kind of failure (format, incorrect data, etc..)
*/
static bool read_data(ifstream &is, grid &b)
{
	unsigned s, i;
	string line;
	bitset<81> coverage;

	while(getline(is, line)) {
		vector<unsigned> g;
		istringstream iss(line);

		if(!(iss >> s)) {
			cerr << "invalid input format" << endl;
			return false;
		}
		while((iss >> i)) {
			if((i < 0) || (i > 80)) {
				cerr << "cell " << i << " out of range" << endl;
				return false;
			}
			if(coverage[i]) {
				cerr << "cell " << i << " multiply covered" << endl;
				return false;
			}

			coverage[i] = true;
			g.push_back(i);
		}
		if(!iss.eof()) {
			cerr << "invalid input format" << endl;
			return false;
		}
		b.add_region(s, g);
	}

	for(i = 0; i < 81; i++)
		if(!coverage[i]) {
			cerr << "uncovered cell " << i << endl;
			return false;
		}
	return true;
}

/*
  Boring but mandatory printing routine. May be upgraded to pretty printing
  some day...
*/
static void print(const grid &b)
{
	for(unsigned i = 0; i < 81; i++) {
		if(i % 9 == 0)
			cout << endl;
		cout << b.cells_[i];
	}
	cout << endl << endl;
}

/*
  Primitive brute-force solver. b is the current state of the grid, and idx
  is the next cell index to fill with a number. Note that the grid state is
  never copied in the search: a reference is always passed. This means that
  placement of a number in the cell must be UNDONE before returning from the
  recursion.
*/
static void solve(grid &b, unsigned idx)
{
	/* if we come this far, the grid contains a solution */
	if(idx > 80) {
		print(b);
		return;
	}

	region *r = b.regions_[idx];
	for(unsigned i = r->lowest_feasible_element(); i < 10; i++) {
		if(b.put(idx, i))
			solve(b, idx+1);
	}

	/* UNDO number placement before returning */
	b.put(idx, 0);
}

int main(int argc, char **argv) try {
	grid b;
	
	if(argc != 2) {
		cerr << "USAGE: " << argv[0] << " INPUT-FILE" << endl;
		return 1;
	}
	ifstream is(argv[1]);
	if(!is) {
		cerr << "can't open input file" << endl;
		return 1;
	}

	if(!read_data(is, b))
		return 1;
	solve(b, 0);
	return 0;
} catch(exception &e) {
	cerr << "INTERNAL ERROR: " << e.what() << endl;
	return 1;
}
