/***************************************************
 * HEADER: Voronoi polyhedron set of programs
 * Vsurfaces, Vcontacts, and Vvolumes, for
 * calculating solvent accessible surfaces,
 * atom-atom contacts, and atom volumes.
 * An analytical solution is used for all
 * calculations.
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CELLSIZE 7.0  // set to (maximum contact radius x 2)
#define PI 3.14159265 // pi.
#define Rw 1.4        // radius of water

// ------------------------- structure definitions ----------------------

struct atom {
  int   atomnum;      // record number from PDB
  float coor[3];      // xyz coordinates of atom
  char  atomname[5];  // name of atom, CA, etc
  char  res[4];       // residue name from PDB
  char chain;
  int   resnum;
  float radius;       // radius of atom
  int   boxnum;       // the box number atom is assigned to.
  char  source;       // source of atom, protein, ligand,
  float SAS;          // solvent exposed surface area
  float vol;          // atom volume
  char  done;         // flag if atom contacts have already been calculated
};

struct atomindex {
  int   nument;      // number of entries in box
  int   first;       // location of first entry in PDBlist
};

struct contactlist {
  int     index;   // index to PDB atom
  double  area;    // contact area, square angstroms
  double  dist;    // distance to atom zero
  char    flag;    // to keep or not. 'X' == omit.
};


// ----------------- Global variables -----------------

struct atom *PDB;           // pointer to PDB array (dynamically allocated)
struct atomindex *box;      // index to PDB atoms within cubic grid
int   *PDBlist;             // list of atoms ordered by box number
struct contactlist contlist[100];   // list of possible atom contacts

// ----------------- prototype definitions ------------------

struct atom *read_PDB(char *FILE_ptr, struct atom *PDB, long int *totatoms);
int    index_protein(int totatoms);
