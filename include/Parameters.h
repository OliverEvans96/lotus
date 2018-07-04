/**
   @file Parameters.h

   @brief Input parameters and parsing functions for the Lotus input file.

   @rst
   See :ref:`input.rst` for a description of the Lotus input file.
   Public member variables of :cpp:class:`Options` are valid keys.
   @endrst
*/


#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"
#include <cassert>
#include <yaml.h>

using namespace std;

typedef map<string, vector<string> > StrVecMap;

///////////////////////////
// High level parameters //
///////////////////////////

/**
   @brief Options from Lotus input file.

   This object is passed to most other objects
   to inform their behavior.
*/
class Options
{
 public:
  /**
     @brief Path to the LAMMPS dumpfile.

     Required.

     @see DumpfileReader
  */
  string dumpfile;

  /**
     @brief Path to the LAMMPS datafile.

     Required.

     @see DatafileReader
  */
  string datafile;

  /**
     @brief Path to output directory.

     Required.

     @see
     @rst
     :ref:`output.rst`
     @endrst
  */
  string outLoc;

  /**
     @brief Number of Timesteps averaged for each frame.

     Required.

     In order to increase the number of atoms in each
     histogram bin, several timesteps are averaged
     together in each frame.
     If the number of steps in the dumpfile is not
     divisible by stepsPerFrame, then the last
     frame will have more steps than the rest.

     @see Time.h
  */
  int stepsPerFrame;

  /**
     @brief List of liquid atom types.

     Required.

     These are the types whose masses contribute towards
     the density calculations for the droplet.

     Lists in YAML can be specified in either of the following ways.
     Both options are equivalent for use in the Lotus input file.

     Option 1:
     ```
     liquidTypes:
       - 1
       - 2
       - 3
     ```

     Option 2:
     ```
     liquidTypes: [1, 2, 3]
     ```

     @see Droplet::fill
     @see Monolayer::fill
  */
  vector<int> liquidTypes;

  /**
     @brief List of susbtrate atom types.

     Required.

     These are the types whose masses contribute towards
     the density calculations for the substrate.
     However, if #substrate is false, feel free to supply
     an empty list (though the option must be present).

     Lists in YAML can be specified in either of the following ways.
     Both options are equivalent for use in the Lotus input file.

     Option 1:
     ```
     substrateTypes:
     - 1
     - 2
     - 3
     ```

     Option 2:
     ```
     substrateTypes: [1, 2, 3]
     ```

     @see Substrate::fill
  */
  vector<int> solidTypes;

  /**
     @brief Number of atoms in the LAMMPS dumpfile.

     Optional. Default: 0

     By default, this is read from the LAMMPS datafile.
     In the case that not all atoms from the datafile
     are present in the dumpfile (e.g. some are frozen),
     it is currently necessary to manually specify
     the number of atoms in the dumpfile via this option.
  */
  int numAtoms;

  /**
     @brief Droplet geometry (spherical or cylindrical).

     Optional. Default: "spherical"
     Acceptable values: "spherical", "cylindrical".

     @rst
     If "spherical" is chosen, then the droplet is assumed
     to be symmetric about the :math:`z` axis, and distance
     to the :math:`z` axis is used as the horizontal axis
     in the 2D density histogram of :ref:`droplet-figure`.

     If "cylindrical" is chosen, then the droplet is assumed
     to be periodic in :math:`y`, and symmetric about the
     :math:`y-z` plane.
     In this case, distance to the :math:`y-z` plane
     is used as the horizontal axis in :ref:`droplet-figure`.
     @endrst

     @see Monolayer::fillOne
     @see Droplet::fillOne
  */
  string geometry; // spherical or cylindrical

  /**
     @brief Whether to save `.png` images of figures.

     Optional. Default: true

     @see
     @rst
     See :ref:`output-image`.
     @endrst
  */
  bool saveImages;

  /**
     @brief Whether to save `.C` ROOT macros of figures.

     Optional. Default: false

     @see
     @rst
     See :ref:`output-root`.
     @endrst
  */
  bool saveROOT;

  /**
     @brief Whether to print verbose output to `stdout`.

     Optional. Default: false

     This may include debugging information and other
     values which are not necessary for routine use.
  */
  bool verbose;

  /**
     @brief LAMMPS bond type for O-H bonds in water molecules.

     Optional. Default: 2

     This is useful if the droplet is water, and
     the orientation of water molecules is calculated.

     @note Molecule orientation is currently not calculated.
     @see DatafileReader::readWaterBonds
  */
  int waterBondType;

  /**
     @brief @f$ z @f$ width of 1D and 2D histogram bins.

     Optional. Default: 1.0

     Units are angstroms.

     @see Grid
  */
  double dz;

  /**
     @brief Volume of 2D histogram bins.

     Optional. Default: 250.0

     Units are cubic angstroms.

     @see Grid
  */
  double dv;

  /**
     @brief Minimum density for TanhFigure and DensFigure.

     Optional. Defaut: 0.0

     Units are g/cm^2

     `TODO`: Should be set to 0.0 everywhere, option is unnecessary.
  */
  double dens_min;

  /**
     @brief Maximum density for TanhFigure and DensFigure.

     Optional. Default: 5.0

     Units are g/cm^3.

     @rst
     `TODO`: Should be combined with :cpp:var:`Options::densMax`.
     The only reason not to combine them is that
     :cpp:class:`DensFigure` (:ref:`dens-figure`) includes
     the substrate, while :cpp:class:`DropletFigure`
     (:ref:`droplet-figure`) does not.
     Since the substrate is generally more dense
     than the liquid, this suggests the use of
     different scales for the two figures.
     @endrst
  */
  double dens_max;

  /**
     @brief
     @rst
     Maximum density for :cpp:class:`DropletFigure` (:ref:`droplet-figure`).
     @endrst

     Optional. Default: 1.5

     Units are g/cm^3.
  */
  double densMax;

  /**
     @brief Maximum @f$ r @f$ value for DropletFigure.

     Optional. Default: 150.0

     This option is not only relevant to the figure produced,
     but specifies the limits of the figure itself and therefore
     the maximum droplet size that can be correctly analyzed.
  */
  double plot_rmax;

  /**
     @brief Maximum @f$ z @f$ value for DropletFigure.

     Optional. Default: 100.0

     This option is not only relevant to the figure produced,
     but specifies the limits of the figure itself and therefore
     the maximum droplet size that can be correctly analyzed.
  */
  double plot_zmax;

  /**
     @brief
     [Aspect ratio](https://en.wikipedia.org/wiki/Aspect_ratio_(image))
     of figures.

     Optional. Default: 1.0

     Plot width is specified by #plot_width, so this option
     is used to implicitly determine plot height.
  */
  double plot_aspect;

  /**
     @brief Width of plots in pixels.

     Optional. Default: 800
  */
  int plot_width;

  /**
     @brief Theoretical density of the bulk liquid droplet.

     Optional. Default: 1.0

     Units are g/cm^3.
     This is useful for estimating the hyperbolic tangent fits
     before actually performing the fits.

     @note This is not actually used, although
     it would be a good idea and easy to implement.
     @see TanhFit::guessTanhFit
  */
  double expectedLiquidDensity;

  /**
     @brief Whether to perform monolayer calculations.

     Optional. Default: true

     By default, the droplet's base radius and contact angle
     are calculated by intersecting the fitted boundary circle
     with the @f$ z @f$ plane defining the top of the monolayer.
     In the case that this option is false, that plane may be
     manually specified through #monoTop.

     @note #monoTop is relative to the top of the substrate.
     @see #substrate
     @see #substrateTop
     @see Monolayer
  */
  bool monolayer;

  /**
     @brief Whether to perform substrate calculations.

     Optional. Default: true

     By default, all @f$ z @f$ coordinates (including #monoTop)
     are given relative to the top of the substrate.
     In the case that this option is false, the reference
     plane may be specified via #substrateTop.

     @see Substrate
  */
  bool substrate;

  /**
     @brief Whether to fit a circle to the boundary points.

     Optional. Default: true

     @rst
     As described in :ref:`overview-quantities`, calculation of
     the bulk radius, contact angle, and bulk height require
     a circle to be fit to the boundary points.
     Therefore, these quantities will not be calculated
     if this option is false.
     @endrst

     @see CircleFit
     @see CircularBulk
  */
  bool fitCircle;

  /**
     @brief @f$ z @f$ coordinate with which to intersect fitted circle.

     Optional. Default: true

     This option is only necessary if #monolayer is false.
  */
  double monoTop;

  /**
     @brief @f$ z @f$ coordinate which is taken to be @f$ z=0 @f$.

     Optional. Default: true

     This option is only necessary if #substrate is false.
  */
  double substrateTop;

  /**
     @brief Radius of cylinder used for 1D droplet density calculations.

     Optional. Default: 10.0

     Units are angstroms.

     In order to calculate a density, we must select a volume
     over which to count atoms and sum their masses.
     For spherical droplets, a thin cylinder around
     the @f$ z @f$ axis is chosen.
     This parameter is the radius of that cylinder.

     @see Droplet::fillOne
     @see Droplet::convertUnits
     @see Droplet::hLiquidDens
  */
  double rDensCyl;

  /**
     @brief Column width of output data files.

     Optional. Default: 15

     Number of characters in each column.
     e.g., the 15 in `%15.6s`.

     @see http://www.cplusplus.com/reference/cstdio/printf/
  */
  int outputColWidth;

  /**
     @brief Floating point precision in output data files.

     Optional. Default: 6

     Number of digits after the decimal point.
     e.g., the 6 in `%15.6s`.

     @see http://www.cplusplus.com/reference/cstdio/printf/
  */
  int outputPrecision;

  /**
     @brief Limits for circle fit parameters.

     Optional.

     Defaults:
     - @f$ x_0 \in [-0.1, 0.1] @f$
     - @f$ y_0 \in [-200.0, 200.0] @f$
     - @f$ r   \in [10.0, 250.0] @f$

     These options limit the potential size and location
     of the fitted circle.

     @see CircleFit::innerFit
  */
  double circleX0Min;
  /// @see #circleX0Min.
  double circleX0Max;
  /// @see #circleX0Min.
  double circleY0Min;
  /// @see #circleX0Min.
  double circleY0Max;
  /// @see #circleX0Min.
  double circleRMin;
  /// @see #circleX0Min.
  double circleRMax;

  void readConfig(const char* configPath);
  void printOptions();
 private:
  StrVecMap yamlMap;
  char configPath[256];

  StrVecMap parseYaml(const char* configPath);
  bool mapHasKey(StrVecMap yamlMap, string key);
  void printYamlMap(StrVecMap yamlMap);
  void fromString(string optionString, bool &option);
  void fromString(string optionString, int &option);
  void fromString(string optionString, double &option);
  void fromString(string optionString, string &option);

  template <typename T>
  void unsafeParseOption(string optionName, T &option);
  template <typename T>
  void unsafeParseOption(string optionName, vector<T> &optionVec);
  template <typename T>
  void parseDefaultOption(string optionName, T &option, T defaultValue);
  template <typename T>
  void parseRequiredOption(string optionName, T &option);

  template <typename T>
  void printOption(string optionName, T option);
  template <typename T>
  void printOption(string optionName, vector<T> option);
};

struct CommandLineParser
{
  CommandLineParser(int argc, const char* argv[]);
  char configPath[256];
  Options options;

  void parseArgs(int argc, const char* argv[]);
  void print();
};

#endif
