/*
  "Kakuro" solver. For details of the game, see e.g.
  http://www.kakuro.com
  
  For instructions on using the program, see the attached README file and
  example inputs.

  Some definitions are in order: RUN: a horizontal or vertical sequence of
  neighbouring cells adding up to a prespecified value. Each empty cell
  belongs to exactly two runs.

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
  In kakuro, the run is always a line (either horizontal or vertical).
  The only pieces of data that we need to describe a run are the number
  of cells and the required sum.

  Since the search procedure iterates only over possible sums for each
  run, there is no need to keep the current sum. If the run is filled,
  it will always add up to the target sum. So we only keep the set of
  numbers currently present in the run,
*/
struct run {
	unsigned	size_;			// number of addends
	unsigned	target_sum_;	// actual sum for reference
	bitset<10> 	set_;			// addends in the set

	run(unsigned target_sum, unsigned size) :
		size_(size), target_sum_(target_sum)
	{ }
	run(unsigned s);

	void add(int n);
};

/*
  This constructor takes the set of included addends as unsigned int
  and based on this calculates the size and the target sum.
*/
run::run(unsigned s) : target_sum_(0), set_(s)
{
	for(s = 1; s < 10; s++)
		if(set_[s])
			target_sum_ += s;
	size_ = set_.count();
}

/*
  Add number n to the set of numbers contained in the run. It should
  never add an existing number. If n is negative, then its value is
  removed from the run.
*/
void run::add(int n)
{
	if(n > 0) {
		if(set_[n])
			throw logic_error("run constraints violated");
		set_[n] = true;
	} else
		set_[-n] = false;
}

/*
  Each run must contain distinct digits 1-9. Therefore:
  - each run has at most 9 cells
  - all sums are composed of individual elements which are in ragne 1-9
  In essence, we have only 2^9=512 distinct sums (some of which may be
  repeated in different places on the grid).

  The table is indexed by the sum value. Each element is a vector of
  run structs; one element for each way the value can be expressed as
  a sum of numbers 1-9.
*/
static vector<vector<run> > G_sum_table;

/* This function populates the G_sum_table. */
static void generate_sum_table()
{
	/* The maximum sum is 45, minimum 0. */
	G_sum_table.resize(46);
	for(unsigned i = 0; i < 512; i++) {
		/* Shift left since 0 is never used in sums. */
		run r(i << 1);
		G_sum_table[r.target_sum_].push_back(r);
	}
}

/* Just a debugging function to print the G_sum_table. */
static void print_sum_table()
{
	unsigned i, j, k;

	for(i = 0; i < 46; i++) {
		cout << "SUM: " << i << endl;
		cout << "==== " << endl;
		for(j = 0; j < G_sum_table[i].size(); j++) {
			cout << G_sum_table[i][j].target_sum_ << " = ";
			for(k = 0; k < 10; k++)
				if(G_sum_table[i][j].set_[k])
					cout << k << " ";
			cout << endl;
		}
		cout << endl;
	}
}

/*
  The grid is rectangular with arbitrary dimensions. Each cell is coded
  according to the input file, according to the following scheme:
  0:		cell to be filled
  -1U:		unfillable cell
  
  Run constraints are not encoded within the grid, but within the run to
  which each cell belongs. Each cell belongs to exactly two runs.
*/
struct grid {
	enum Direction {
		DOWN = 0, RIGHT
	};

	struct run_pair {
		run *r_[2];

		run_pair(run *d = 0, run *r = 0)
		{
			r_[0] = d; r_[1] = r;
		}
	};
	
	unsigned		rows_, cols_;	// grid dimensions
	vector<unsigned> cells_;		// all cells
	vector<run_pair> runs_;			// cell -> run map

	grid() : rows_(0), cols_(0)
	{ }

	void set_size(unsigned rows, unsigned cols);
	bool put(unsigned idx, unsigned n);
	bool add_run(unsigned startidx, unsigned target_sum, unsigned size, Direction dir);
};

/*
  Set size of the grid and initialize cells to appropriate values so that we
  can check for validity of input data.

  Cell values are initialized to -2 so that we can detect uncovered cells,
  unfillable cells and empty cells.
*/
void grid::set_size(unsigned rows, unsigned cols)
{
	rows_ = rows; cols_ = cols;

	cells_.resize(rows_ * cols_, -2U);
	runs_.resize(rows_ * cols_);
}

/*
  Put number n into cell indexed by idx. If no constraints are violated, the
  method puts the number in the cell (updating its corresponding run) and
  returns true, otherwise it returns false. Putting 0 clears the cell,
  whatever its contents.
*/
bool grid::put(unsigned idx, unsigned n)
{
	unsigned oldn = cells_[idx];
	
	if(n > 9)
		throw invalid_argument("invalid argument to grid::put()");

	/* check that the number doesn't exist in any run */
	if(runs_[idx].r_[0]->set_[n] || runs_[idx].r_[1]->set_[n])
		return false;

	/* remove any previous number */
	runs_[idx].r_[0]->add(-oldn);
	runs_[idx].r_[1]->add(-oldn);

	/* add new number */
	runs_[idx].r_[0]->add(n);
	runs_[idx].r_[1]->add(n);
	cells_[idx] = n;

	return true;
}

/*
  Add a new run to the grid. The run is described with the following:
  - starting position (the index of the CONSTRAINT cell; it is also added
    to unfillable cells)
  - required sum
  - size (number of cells)
  - dir: down (0) or right (1)

  Returns false if we try to cover an already covered cell.
*/
bool grid::add_run(
	unsigned startidx,
	unsigned target_sum,
	unsigned size,
	Direction dir)
{
	unsigned incr = dir == RIGHT ? 1 : cols_;

	if((dir != DOWN) && (dir != RIGHT))
		throw invalid_argument("invalid arguments to grid::add_run()");

	/*
	  TODO! check here for startidx overflow before indexing!!

	  constraint cell must fall into uninitialized cells or onto another
	  constraint cell (giving the other part of the constraint). so we
	  accept all values (-1 and -2) except 0.
	*/
	if(!cells_[startidx])
		return false;

	/* add constraint cell to unfillable cells */
	cells_[startidx] = -1U;

	/* create run and update back-references */
	run *r = new run(target_sum, size);
	for(startidx += incr; size; startidx += incr, size--) {
		/*	
		  check for grid limits. this is not comprehensive check as
		  the specified run might run over the right edge of the
		  grid. but it should be caught by the latter logic as the
		  left grid is unfillable.
		*/
		if(startidx >= rows_ * cols_)
			return false;
	
		/* check for overflow into unfillable space */
		if(cells_[startidx] == -1U) 
			return false;
		
		/* check for cell occupied by another run in the same direction */
		if(runs_[startidx].r_[dir])
			return false;

		/* assign the cell to a run and mark it as empty */
		runs_[startidx].r_[dir] = r;
		cells_[startidx] = 0;
	}

	return true;
}

/*
  Load data, performing sanity checking along the way. The data format is
  described in the README file. The method returns true on success, and false
  on any kind of failure (format, incorrect data, etc..)
*/
static bool read_data(ifstream &is, grid &b)
{
	unsigned idx;
	string line;

	if(!getline(is, line)) {
		cerr << "error reading the dimensions line" << endl;
		return false;
	} else {
		unsigned rows, cols;
		istringstream iss(line);

		if(!(iss >> rows >> cols)) {
			cerr << "invalid dimensions line" << endl;
			return false;
		}

		/* it doesn't make sense to have less that 2x2 puzzle */
		if((rows < 2) || (cols < 2)) {
			cerr << "invalid dimensions (must be at least 2 in each direction)" << endl;
			return false;
		}

		b.set_size(rows, cols);
	}

	/* process unfillable cells */
	if(!getline(is, line)) {
		cerr << "error reading the unfillable cells line" << endl;
		return false;
	} else {
		istringstream iss(line);
		
		while(iss >> idx) {
			// TODO: range check
			b.cells_[idx] = -1U;
		}
		if(!iss.eof()) {
			cerr << "invalid unfillable cells line" << endl;
			return false;
		}
	}

	/* process constraint lines */
	while(getline(is, line)) {
		unsigned sum[2], size[2];
		istringstream iss(line);

		if(!(iss >> idx >> sum[0] >> size[0] >> sum[1] >> size[1])) {
			cerr << "invalid input format for constraint index " << idx
				 << endl;
			return false;
		}
		if((sum[0] > 45) || (sum[1] > 45) || (size[0] > 9) || (size[1] > 9)) {
			cerr << "inconsistent data for constraint index " << idx << endl;
			return false;
		}
		if(!b.add_run(idx, sum[0], size[0], grid::DOWN)) {
			cerr << "overflowing (down?) run or misplaced constraint cell " << idx << endl;
			return false;
		}
		if(!b.add_run(idx, sum[1], size[1], grid::RIGHT)) {
			cerr << "overflowing (right?) run or misplaced constraint cell " << idx << endl;
			return false;
		}
	}

	/* check that each 0-cell is a member of exactly two runs */
	for(idx = 0; idx < b.cells_.size(); idx++) {
		if(b.cells_[idx] == -2U) {
			cerr << "uncovered cell " << idx << endl;
			return false;
		}

		if((b.cells_[idx] != 0) && (b.cells_[idx] != -1U))
			throw logic_error("invalid grid initialization");

		if((b.cells_[idx] == 0)
		&& (!b.runs_[idx].r_[0] || !b.runs_[idx].r_[1])) {
			cerr << "incompletely covered cell " << idx << endl;
			return false;
		}
	}

	return true;
}

/*
  Boring but mandatory printing routine. May be upgraded to pretty printing
  some day...
*/
static void print(const grid &b)
{
	for(unsigned i = 0; i < b.rows_ * b.cols_; i++) {
		if(i % b.cols_ == 0)
			cout << endl;
		if(b.cells_[i] == -1U)
			cout << '.';
		else
			cout << b.cells_[i];
	}
	cout << endl << endl;
}

/*
  Return a set of feasible addends for the given run. It doesn't return
  numbers that are already present in the run.

  In order to generate correct sums, we choose only those possible sums
  which contain ALL elements already present in the run.
*/
static bitset<10> feasible_addend_set(const run *r)
{
	bitset<10> result;
	
	if(r->target_sum_ > 45)
		throw invalid_argument("run inconsistency detected in feasible_addend_set()");

	/*
	  TODO: This loop can be avoided by secondary index on the number
	  of addends in the G_sum_table.
	*/
	for(unsigned i = 0; i < G_sum_table[r->target_sum_].size(); i++) {
		run ns = G_sum_table[r->target_sum_][i];
		
		/*
		  the chosen sum must have the same number of addends as the size of
		  the run, and the current contents of the run must be a SUBSET of
		  the chosen sum. the latter condition guarantees that the sum in
		  the filled run will equal the target sum.
		*/
		if((ns.size_ == r->size_) && ((ns.set_ & r->set_) == r->set_))
			result |= ns.set_;
	}

	/* remove from the feasible set numbers already present in the run */
	return result & ~r->set_;
}

/*
  Find the set of feasible numbers to put in the cell at index idx.

  TODO: can be made faster by removing if()s and pointing to empty sets
  within the grid itself.
*/
static bitset<10> feasible_numbers(const grid &b, unsigned idx)
{
	bitset<10> result_down, result_right;

	/* Find the feasible numbers for the down and right runs. */
	if(b.runs_[idx].r_[0])
		result_down = feasible_addend_set(b.runs_[idx].r_[0]);
	if(b.runs_[idx].r_[1])
		result_right = feasible_addend_set(b.runs_[idx].r_[1]);

	/* Since the cell is member of both runs, INTERESECT feasible sets. */
	return result_down & result_right;
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
	if(idx >= b.rows_ * b.cols_) {
		print(b);
		return;
	}

	/*
	  for coding simplicity, we always recurse, even on nonfillable cells.
	  so skip them here..
	*/
	if(b.cells_[idx] == -1U) {
		solve(b, idx+1);
		return;
	}

	bitset<10> trial_numbers = feasible_numbers(b, idx);
	for(unsigned i = 1; i < 10; i++)
		if(trial_numbers[i]) {
			/*
			  putting should always succeed since the trial number set is
			  calculated so that it contains only not present numbers.
			*/
			if(!b.put(idx, i))
				throw logic_error("duplicate trial number");
			solve(b, idx+1);
		}

	/* UNDO number placement before backtracking */
	b.put(idx, 0);
}

int main(int argc, char **argv) try {
	grid b;
	
	if(argc != 2) {
		cerr << "USAGE: " << argv[0] << " INPUT-FILE" << endl;
		return 1;
	}

	generate_sum_table();
	if(string(argv[1]) == "-test") {
		print_sum_table();
		return 0;
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
