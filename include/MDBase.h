/**
   @file MDBase.h

   @brief Basic quantities related to the MD simulation.

   Contains time-independent data such as grid and
   atom information, as well as time-dependent
   information such as the substrate and monolayer extents.
*/

#ifndef MDBASE_H
#define MDBASE_H

#include "Utils.h"
#include "Parameters.h"
#include "Time.h"

using namespace std;

/////////////////////
// Simulation Data //
/////////////////////

/**
  @brief Simulation box bounds, as defined in the LAMMPS dumpfile.

  This data is updated each timestep.
  @see DumpfileReader
*/
struct BoxBounds {
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;
};

/**
  @brief Source of truth for general MD variables.

  Variables are set here, and a pointer to this object
  is passed around and read elsewhere.
*/
struct SimData {
  // @brief Static data
  Options options;
  int numAtoms;
  int numSteps;
  int numFrames;
  int stepsPerFrame;
  LastFrame lastFrame;
  /// @brief Which atom types correspond to liquid
  vector<int> liquidTypes;
  /// @brief Which atom types correspond to substrate
  vector<int> solidTypes;
  /// @brief Mass for each atom type
  map<int, double> masses;
  /** @brief Maps oxygen atom id to two hydrogen atom ids.
      @see DatafileReader::readWaterBonds
  */
  map<int, int*> waterBonds;
  BoxBounds simBounds;

  // Dynamic data
  /** @brief z coordinate of the top of the substrate.

     If the substrate option is disabled,
     this can be manually specified via Options::substrateTop.
     @see Substrate::findLimits
  */
  double substrateTop;
  /** @brief z coordinate of the top of the monolayer.

     If the monolayer option is disabled,
     this can be manually specified in Options::monoTop.
     @see Droplet::findMonolayer
     @see Monolayer::findMonoLimits
  */
  double monoTop;
  Frame *framePtr;

  SimData(Options options);
  ~SimData();

  /** @brief Since water bonds are stored as a map from int to int*,
      those int*s must be allocated upon storage in the map.
      So we delete them here.

      @see DatafileReader::readWaterBonds
  */
  void deleteWaterBonds();
  void setOptions(Options options);
  void setNumSteps(int _numSteps);
  void setStepsPerFrame(int _stepsPerFrame);
};

/**
   @brief Grid to use for 2D histogram.

   Each bin in the grid has the same volume.
   #dz and #dv are specified by Options::dz and Options::dv,
   and determine the size of the bins.
   `dr` is determined automatically for each bin from the other two.
*/
struct Grid {
  double zlo, zhi, rhi, vhi;
  double dz, dv;
  int nz, nv, nr;

  /** @brief @f$r @f$ values of bin edges (including both min and max)
      @see #calculateBins
      @see #allocateBins
  */
  double* rVals;
  /// @see calculateVolumeLimits
  double zhi_preround;
  /// @see calculateVolumeLimits
  double vhi_preround;
  /// @see calculateVolumeLimits
  double rhi_preround;
  /// @brief Whether #rVals has been allocated
  bool allocated;

  ~Grid();

  void setBounds(double _zlo, double _zhi, double _rhi);
  void setSpacing(double _dz, double _dv);

  /// @brief Set up grid after calling #setBounds and #setSpacing
  void init();
  /** @brief Round upper z and r bounds in case #dv or #dz
      don't evenly divide the z and r ranges.
  */
  void calculateVolumeLimits();
  /// @brief Allocate #rVals before calling #calculateBins
  void allocateBins();
  /// @brief Set #rVals after calling #allocateBins
  void calculateBins();
  /// @brief Deallocate #rVals after using
  void deallocateBins();
};

#endif
