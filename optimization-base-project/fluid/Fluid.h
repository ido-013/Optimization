
#pragma once

#pragma warning(push, 0)
#include <d3dx9math.h>
#pragma warning(pop)
#include <list>
#include <vector>

// Fluid magic numbers
const float FluidTimestep = 0.005f;
const float FluidSmoothLen = 0.012f;
const float FluidStaticStiff = 3000.0f;
const float FluidRestDensity = 1000.0f;
const float FluidWaterMass = 0.0002f;
const float FluidViscosity = 0.1f;
const float FluidStiff = 200.0f;
const float FluidInitialSpacing = 0.0045f;

/*****************************************************************************/

class Fluid 
{
	public:
		/* Common Interface */
		Fluid();
		~Fluid();

		void Create(double w, double h);
		void Fill(double size);
		void Clear();
		void Update(double dt);

		/* Common Data */
		unsigned int * gridindices;

		std::vector<D3DXVECTOR2> particle_positions;
		std::vector<D3DXVECTOR2> particle_velocities;
		std::vector<D3DXVECTOR2> particle_accelerations;
		std::vector<double> particle_densities;
		std::vector<double> particle_pressures;

		std::vector<unsigned int> gridoffsets_offset;
		std::vector<unsigned int> gridoffsets_count;

		unsigned int neighbors_capacity;
		unsigned int num_neighbors;

		std::vector<unsigned int> neighbors_p;
		std::vector<unsigned int> neighbors_n;
		std::vector<double> neighbors_distsq;

		unsigned int Size()					{ return particles_size; }
		unsigned int Step()					{ return step; }
		void Pause( bool p )				{ paused = p; }
		void PauseOnStep( unsigned int p )	{ pause_step = p; }
		double Width()						{ return width; }
		double Height()						{ return height; }

	private:
		
		/* Simulation */
		void UpdateGrid();
		void GetNeighbors();
		void ComputeDensity();
		void SqrtDist();
		void ComputeForce();
		void Integrate(double dt);

	private:
		/* Run State */
		unsigned int step;
		bool paused;
		unsigned int pause_step;

		/* World Size */
		double width;
		double height;
		int grid_w;
		int grid_h;

		unsigned int particles_size;

		/* Coefficients for kernel */
		double poly6_coef;
		double grad_spiky_coef;
		double lap_vis_coef;
};
