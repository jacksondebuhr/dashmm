#ifndef __DASHMM_DISTRIBUTER_H__
#define __DASHMM_DISTRIBUTER_H__
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <float.h>
#include <math.h>
#include <assert.h> 

// Cell holds information on how it is partitioned and how it relates to 
// its spacial neighbors and the other cells in its partition 
struct Cell{
	int x;
	int y;	
	int z;
	int nparticles;
	int available_neighbors;
	int group_neighbors;
	double structure_priority;
	double neighborness;
	double priority;
	bool is_free;
	int part_index;

	Cell() : part_index(-1) is_free(true) {}
};

class Distributer{
public:
	Distributer() {}

	// 
	void setup(int nparts, int*** counts, int len) {
		length = (int)std::pow(len, 1/3);
		nparticles = 0;
		npartitions = nparts;
		std::vector<std::vector<std::vector<Cell*> > > g;
		for(int x=0;x<length;++x){
			std::vector<std::vector<Cell*> > base;
			for(int y=0;y<length;++y){
				std::vector<Cell*> col;
				for(int z=0;z<length;++z){
					Cell* cell = new Cell;
					cell->x = x; 
					cell->y = y;
					cell->z = z;
					cell->num_particles = counts[x][y][z];
					num_particles += counts[x][y][z];
					col.push_back(cell);	
				}
				base.push_back(col);
			}
			g.push_back(base);
		}
		grid = g;
	}

	int* partition()


private:
	// Partitioning Algorithms


	// REPRESENTATION
	int length; 
	int nparticles;
	int npartitions;
	std::vector<std::vector<std::vector<Cell*> > > grid;
	std::vector<std::vector<Cell*> > partitions;

};




















#endif
