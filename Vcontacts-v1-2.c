
/**********************************************************************
 * PROGRAM Vcontacts.c
 * calculates the atom-atom contacts within a protein, along with the
 * solvent accessible surface.
 *
 * Please cite the following article as a reference:
 * Quantification of protein surfaces, volumes and atom-atom contacts using a constrained Voronoi procedure
 * B.J.McConkey, V. Sobolev, and M. Edelman (2002), Bioinformatics, 18:1365-1373.
 * ====================================================================

 * NOTES ON METHODOLOGY:
 *
 * The method calculates the planes of contact the center atom makes 
 * will each atom in the contact list, then determines the points
 * of intersection of the planes and intersections with a sphere. The
 * sphere has a radius equal to the sum of the van der Waals radius of the 
 * center atom plus the radius of a solvent atom (water). The set of 
 * intersections of the planes defines a polygon surrounding the center 
 * atom. This polygon is projected onto the surface of the sphere. The
 * area of each projection is calculated as the sum of spherical triangles
 * and arc segments, where an arc segment is the difference between a 
 * angluar segment of a spherical cap and the corresponding spherical triangle.
 *
 * The points of intersection are calculated using a convex hull algorithm.
 * The algorithm operates with efficiency O(Nk), where N is the number of 
 * input atoms and k is the number of edges on the contact polyhedron.
 *
 * The convex hull algorithm effectively describes the atom contacts only if
 * no 'engulfing' atoms are present, ie., atoms which cover more than 50% of
 * the atom contact surface. This may occur between closely bonded atoms of
 * different radii, such as some C=O bonds. A correction factor accounts for 
 * this.
 **********************************************************************/

#include "Vpoly.h"

#define PI 3.14159265       // pi.
#define Rw 1.4              // radius of water
#define PAUSE ch=getchar()  // for debugging

// ------------------ structure definitions --------------------

struct plane {
  double Ai[4];      // parameters A,B,C,D of contact plane Ax+By+Cz+D=0 
  double dist;       // distance from plane to origin
  int    index;      // index to which record in PDB or ligand array.
  double area;       // contact area in square angstroms
  char   flag;       // 'X' if no contact, 'E' if an engulfing atom.
};

struct vertex {
  double xi[3];      // x,y,z coordinates (x1,x2,x3)
  double dist;       // distance to origin
  int    plane[3];   // identification of intersecting planes. -1 = sphere.
};

struct ptindex {
  int  numpts;       // number of points defining face
  int  pt[40];       // index to polyhedron points
};

struct edgevector {
  double V[3];       // vector for edge
  int startpt;       // initial vertex
  int endpt;         // final vertex
  int plane[2];      // planes defining edge (using single atoms for now)
  int startplane;    // third plane at start point
  int endplane;      // third plane at end point
  char arc;          // flag for arc point calculations
};

struct ca_struct {
  long int prev;     // previous contact location in ca_index
  long int atom;     // PDBarray number of current contact (NOT PDB record number)
  float area;        // contact area
  float dist;        // distance between atoms
};

// ----------------- global variables -------------------

struct ptindex    ptorder[100];  // for ordering vertices around each face
struct vertex     centerpt[100]; // center points for each contact face
struct vertex     poly[200];     // polyhedron vertices
struct plane      cont[100];     // atom and contact plane information
struct edgevector vedge[200];
struct ca_struct *ca_rec;        // array - contact area records
long int *ca_index;              // array - index to first ca_recs for each atom.
long int  numcarec = 0;          
long int  ca_recsize;
int       dim;                   // size of side of box, in units CELLSIZE.
float     globalmin[3];          // minimum coordinate value for protein atoms
float     globalmax[3];          // maximum coordinate value for protein atoms
char      ch;                    // for debugging 
char      planedef;              // = X, R, or B
char      showbonded;            // = Y or N.
char      normalize;             // = Y or N. normalize areas to area of sphere.
long int *seed;                  // seed vertices for new polyhedra

// --------------------- function prototypes ----------------------

double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado);
double cosPQR( double ptP[], double ptQ[], double ptR[]);
char   test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC);
char   order_faces(int atomzero, struct vertex poly[], struct vertex centerpt[], float rado, int NC, 
		   int NV, struct plane cont[], struct ptindex ptorder[]);
void   project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC, 
		      int NV, struct plane cont[]);
int    voronoi_poly2(struct atom PDB[], int atomzero, struct plane cont[], float rado, 
		     int NC, struct contactlist contlist[]);
int    add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC);
void   add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1, 
		 int testplane, struct vertex poly[], int startpt);
int    solve_3x3(double eq0[], double eq1[], double eq2[], double pt[]);
int    solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[]);
int    calc_region(struct atom PDB[], int PDBtot, int atomlist[], int numatoms);
void   calc_areas(struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, 
		  struct plane cont[], struct ptindex ptorder[], int atomzero);
int    index_protein(int PDBtot);
void   assign_radii(struct atom *PDB, long int PDBtot);
void   save_areas(struct plane cont[], struct contactlist contlist[], int NC, int atomzero);
int    get_contlist4(struct atom PDB[], int atomzero, struct contactlist contlist[], 
		     int PDBtot, float rado, int dim);
void   parse_commandline(int argc, char *argv[], char FILENAME[]);
void   save_seeds(struct plane cont[], struct vertex poly[], int NV, int atomzero);
void   get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC, int atomzero);


// ========================================================================

int main(int argc, char *argv[])
{
  int      atomi;           // atom number of atomzero
  char     FILENAME[500];    // input file name
  long int PDBtot;          // number of PDB atoms

  parse_commandline( argc, argv, FILENAME);
  PDB = read_PDB(FILENAME, PDB, &PDBtot);

  // initialize contact atom index
  ca_index = malloc(PDBtot*sizeof(long int));
  ca_recsize = 5*PDBtot;
  ca_rec = malloc(ca_recsize*sizeof(struct ca_struct));
  seed = malloc(3*PDBtot*sizeof(long int));
  if((!ca_rec) || (!ca_index) || (!seed)) {
    printf("memory allocation error\n"); 
    exit(1);
  }
  for(atomi=0; atomi<PDBtot; ++atomi) {
    PDB[atomi].vol = 0.0;
    ca_index[atomi] = -1;   //initialize pointer array
    seed[atomi*3] = -1;     // initialize seed array
  }

  // assign protein atoms to boxes in cubic grid
  dim = index_protein(PDBtot);

  // calc volumes for all protein atoms
  calc_region(PDB, PDBtot, PDBlist, PDBtot);
  return(0);
}

/***************************
 * subroutine calc_region
 ***************************/

// this subroutine calculates the contact areas (SAS) for a given set of atoms.
// Here, the set of atoms is the entire protein.
// needs global variable 'dim'.

int calc_region(struct atom PDB[], int PDBtot, int atomlist[], int numatoms)
{
  int    atomi;    // atom counter
  int    atomzero; // current center atom
  int    cai;      // contact atom counter for atomzero
  int    NC;       // number of contacts around atomzero
  int    NV;       // number of vertices in polyhedron around atomzero
  double areao;    // area of atom zero
  double areatot;  // total contact area
  float  rado;     // radius of atomzero PLUS radius of water
  float  SAS;      // solvent exposed surface in square angstroms
  char   surfatom; // atom type, 'I' internal, 'S' surface
  float  SAStot;
  long int currindex;
  int      contnum;
  float *coorA, *coorB; // coordinate pointers

  for(atomi=0; atomi<numatoms; ++atomi) {

    // ============= atom contact calculations =============
    atomzero = atomlist[atomi];
    rado = PDB[atomzero].radius + Rw;
    NC = get_contlist4(PDB, atomzero, contlist, PDBtot, rado, dim);
    NV = voronoi_poly2(PDB, atomzero, cont, rado, NC, contlist);
    surfatom = order_faces(atomzero, poly, centerpt, rado, NC, NV, cont, ptorder);
    calc_areas(poly, centerpt, rado, NC, NV, cont, ptorder, atomzero);
    save_areas(cont, contlist, NC, atomzero);

    // ============  SAS contact area calculations ============
    areao = 4.0*PI*rado*rado;
    areatot = 0.0;
    for(cai=0; cai<NC; ++cai) {
      if(cont[cai].flag != 'X') {
	areatot += cont[cai].area;
      }
    }
    SAS = areao - areatot;
    if( SAS > 0.0) {  // remove fractional negative values (~ -0.0000002)
      PDB[atomzero].SAS = SAS;
      SAStot += SAS;
    } else {
      PDB[atomzero].SAS = 0.0;
    }
  }

  // ==================== write output file =======================

  // -------------- print header ---------------
  printf("#  Atom-atom contact areas, calculated using program 'Vcontacts'\n");
  printf("#  version 1.2\n");
  printf("#  Authors: B.J.McConkey, V.Sobolev, and M.Edelman\n");
  printf("#  \n");
  if(showbonded == 'Y') {
    printf("#    -all: contacts include covalently bonded atoms\n");
  } else {
    printf("#    covalently bonded atom contacts are not listed\n");
  }
  if(planedef == 'X') {
    printf("#    -planedef X: extended radical plane\n");
  } else if(planedef == 'R') {
    printf("#    -planedef R: radical plane\n"); // used as default
  } else if(planedef == 'B') {
    printf("#    -planedef B: bisecting dividing plane\n");
  }
  if(normalize == 'Y') {
    printf("#    -norm: contacts normalized to a percent of total contact area\n");
  } else {
    printf("#    Contacts are given in SAS equivalent units (square angstroms)\n");
  }
  printf("#   \n");
  printf("#  -------ATOM-------   ------CONTACT----- \n");
  printf("#  NUM NAME   RESIDUE   NUM NAME   RESIDUE   AREA%c  DIST \n", 
	 (normalize == 'Y')?'%':' ');
  printf("#--------------------------------------------------------\n");


  // print SAS and contact areas
  for(atomi=0; atomi<PDBtot; ++atomi) {
    currindex = ca_index[atomi];
    rado = PDB[atomi].radius+Rw;
    contnum = 1;
    printf("%5d %5s %4d %5s %1c       Sol_acc_surf      %7.2f\n", PDB[atomi].atomnum, 
	   PDB[atomi].atomname, PDB[atomi].resnum, PDB[atomi].res, PDB[atomi].chain,
	   (normalize=='Y')?(PDB[atomi].SAS*100.0/(rado*rado*4.0*PI)):PDB[atomi].SAS);
    while(currindex != -1) {
      // normalize if required
      if(normalize == 'Y') {
	ca_rec[currindex].area *= (100.0/(rado*rado*4.0*PI));
      }
      if(showbonded == 'Y') {
	printf("                      %5d %5s %4d %5s %1c %7.2f  %4.2f\n", 
	       PDB[ca_rec[currindex].atom].atomnum, PDB[ca_rec[currindex].atom].atomname, 
	       PDB[ca_rec[currindex].atom].resnum, PDB[ca_rec[currindex].atom].res, PDB[atomi].chain, 
	       ca_rec[currindex].area, ca_rec[currindex].dist);

	++contnum;
      } else {
	// use distance check d<2.0 A to identify bonded contacts. always print S-S bonds.
	coorA = PDB[atomi].coor;
	coorB = PDB[ca_rec[currindex].atom].coor;
	if((((coorA[0]-coorB[0])*(coorA[0]-coorB[0]) + (coorA[1]-coorB[1])*(coorA[1]-coorB[1]) 
	    + (coorA[2]-coorB[2])*(coorA[2]-coorB[2])) > 4.0) 
	   || ((PDB[atomi].atomname[1] == 'S') && (PDB[ca_rec[currindex].atom].atomname[1] == 'S')))  {
	  printf("                      %5d %5s %4d %5s %1c %7.2f  %4.2f\n", 
		 PDB[ca_rec[currindex].atom].atomnum, PDB[ca_rec[currindex].atom].atomname, 
		 PDB[ca_rec[currindex].atom].resnum, PDB[ca_rec[currindex].atom].res, PDB[atomi].chain, 
		 ca_rec[currindex].area, ca_rec[currindex].dist);
	  ++contnum;
	}
      }
      currindex = ca_rec[currindex].prev;
    }
    printf("\n");
  }
  return(0);
}


/******************************
 * subroutine voronoi_poly2
 * created 08/07/2001  BJM
 ******************************/

int voronoi_poly2(struct atom PDB[], int atomzero, struct plane cont[], float rado, 
		  int NC, struct contactlist contlist[])
{
  int    cai;          // contact atom counter
  struct atom *ca_ptr; // pointer to pdb atom
  double atomdist;     // distance to atom
  double planedist;    // distance to plane
  double mindist;      // distance to closest plane
  int    planeA;       // closest plane to origin
  int    planeB;       // second plane, with planeA defines closest edge
  int    planeC;       // new intersection plane for edge (endpt)
  int    oldplaneC;    // old intersection plane for edge (startpt)
  double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
  double vtmin;        // minimum value for vt 
  double vtdiv;        // check for division by zero in vt calculation 
  double temppt[3];
  struct vertex *stp; // pointer to start vertex, coordinates
  double *V;          // pointer to edge vector
  int    startedge;
  int    edgenum;
  int    vn = 0;
  char   edgeflag;
  int    edgei;      // edge counter
  int    vi, vj;     // vertices counters
  double arcpt0[3], arcpt1[3];
  int    testpA, testpB;
  double testvalA, testvalB;
  char   arcflag = 'N';

  // failsafe variables:
  char   recalc;       // flag if hull is being recalculated (orig. unbounded)
  float  origcoor[3];  // original pdb coordinates for atom. 

  recalc = 'N';
 RESTART:
  planeA = -1;
  planeB = -1;
  planeC = -1;

  /* generate planes of contact with D = planedist */
  mindist = 9.9e+9;
  for(cai=0; cai<NC; ++cai) {
    ca_ptr = &PDB[contlist[cai].index];
    atomdist = contlist[cai].dist;
    
    if(planedef == 'B') {  // bisection - original Voronoi procedure
      planedist = atomdist/2.0;
    } else if(planedef == 'R') { // radical plane (Gellatly and Finney) - default.
      planedist = (atomdist*atomdist + (rado-Rw)*(rado-Rw) - (ca_ptr->radius)*(ca_ptr->radius))/(2*atomdist);
    } else { // extended radical plane (McConkey et al). 
      planedist = (atomdist*atomdist + rado*rado - (Rw + ca_ptr->radius)*(Rw + ca_ptr->radius))/(2*atomdist);
    } 

    cont[cai].Ai[0] = (ca_ptr->coor[0] - PDB[atomzero].coor[0])/atomdist;
    cont[cai].Ai[1] = (ca_ptr->coor[1] - PDB[atomzero].coor[1])/atomdist;
    cont[cai].Ai[2] = (ca_ptr->coor[2] - PDB[atomzero].coor[2])/atomdist;
    cont[cai].Ai[3] = -planedist;
    cont[cai].dist  = fabs(planedist);
    cont[cai].index = contlist[cai].index;
    cont[cai].flag = 'X'; // initialize contact flags to 'no contact' 

    // set plane0 as closest plane
    if(cont[cai].dist < mindist) {
      mindist = cont[cai].dist;
      planeA = cai;
    }
  }

  // add four planes surrounding atom, outer limit for voronoi polyhedron
  cont[NC].Ai[0] = 0.707;
  cont[NC].Ai[1] = 1.0;
  cont[NC].Ai[2] = 0.0;
  cont[NC].Ai[3] = -10.0;
  cont[NC+1].Ai[0] = 0.707;
  cont[NC+1].Ai[1] = -1.0;
  cont[NC+1].Ai[2] = 0.0;
  cont[NC+1].Ai[3] = -10.0;
  cont[NC+2].Ai[0] = -0.707;
  cont[NC+2].Ai[1] = 0.0;
  cont[NC+2].Ai[2] = 1.0;
  cont[NC+2].Ai[3] = -10.0;
  cont[NC+3].Ai[0] = -0.707;
  cont[NC+3].Ai[1] = 0.0;
  cont[NC+3].Ai[2] = -1.0;
  cont[NC+3].Ai[3] = -10.0;

  // get starting vertex from seed or calc new vertex
  get_firstvert(cont, &planeA, &planeB, &planeC, NC, atomzero);
  solve_3x3(cont[planeA].Ai, cont[planeB].Ai, cont[planeC].Ai, temppt);

  // add first vertex to vertex list
  add_vertex(poly, 0, temppt, planeA, planeB, planeC);
  
  // flag contacts as present
  cont[planeA].flag = 'Y';
  cont[planeB].flag = 'Y';
  cont[planeC].flag = 'Y';
  
  // calculate edge vectors
  add_vedge(vedge, 0, cont, planeA, planeB, planeC, poly, 0);
  add_vedge(vedge, 1, cont, planeB, planeC, planeA, poly, 0);
  add_vedge(vedge, 2, cont, planeC, planeA, planeB, poly, 0);
  
  startedge = 0;
  edgenum = 3;
  vn = 1;

  /* --------------------------------------------------- */
  /* Generate new polyhedron points from edge vectors    */
  /* --------------------------------------------------- */

  while(1) {
    // get next unfinished vector = startedge
    while((vedge[startedge].endpt >= 0) && ((edgenum-startedge) > 0)) {
      ++startedge;
    }
    if((edgenum-startedge) <= 0) {
      // all edges are done, polyhedron complete.
      break;
    }

    vtmin = 9.9e+9; // dummy value
    stp = &poly[vedge[startedge].startpt];
    V = vedge[startedge].V;
    planeA = vedge[startedge].plane[0];
    planeB = vedge[startedge].plane[1];
    oldplaneC = vedge[startedge].startplane;
    planeC = -1;

    // get closest positive intersection point
    for(cai=0; cai<NC; ++cai) {
      // check if contact is to be done - for now, do all.
      if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
	vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
	if(vtdiv != 0.0) {
	  vt = -(cont[cai].Ai[0]*stp->xi[0] +cont[cai].Ai[1]*stp->xi[1] 
		 +cont[cai].Ai[2]*stp->xi[2] +cont[cai].Ai[3])/vtdiv;
	  if((vt < vtmin) && (vt > 0)) {
	    vtmin = vt;
	    planeC = cai;
	  }
	}
      }
    }
    poly[vn].xi[0] = stp->xi[0] + vtmin*V[0];
    poly[vn].xi[1] = stp->xi[1] + vtmin*V[1];
    poly[vn].xi[2] = stp->xi[2] + vtmin*V[2];

    // if point is outside sphere, check vs. external planes
    if((poly[vn].xi[0]*poly[vn].xi[0] + poly[vn].xi[1]*poly[vn].xi[1] 
	+ poly[vn].xi[2]*poly[vn].xi[2]) > rado*rado) {
      for(cai=NC; cai<NC+4; ++cai) {
	// check if contact is to be done - for now, do all.
	if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
	  vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
	  if(vtdiv != 0.0) {
	    vt = -(cont[cai].Ai[0]*stp->xi[0] +cont[cai].Ai[1]*stp->xi[1] 
		   +cont[cai].Ai[2]*stp->xi[2] +cont[cai].Ai[3])/vtdiv;
	    if((vt < vtmin) && (vt > 0)) {
	      vtmin = vt;
	      planeC = cai;
	    }
	  }
	}
      }
      poly[vn].xi[0] = stp->xi[0] + vtmin*V[0];
      poly[vn].xi[1] = stp->xi[1] + vtmin*V[1];
      poly[vn].xi[2] = stp->xi[2] + vtmin*V[2];
    }

    add_vertex(poly, vn, poly[vn].xi, planeA, planeB, planeC);
    vedge[startedge].endpt = vn;
    vedge[startedge].endplane = planeC;

    //flag contact as present
    cont[planeC].flag = 'Y';

    // ========  ADD EDGES  ========

    // check edge (planeA, planeC)
    edgeflag = 'Y';
    edgei = startedge+1;
    while(edgei < edgenum) {
      if(((vedge[edgei].plane[0] == planeA)&&(vedge[edgei].plane[1] == planeC)) ||
	 ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeA))) {
	// already on list, add current vertex as endpt
	vedge[edgei].endpt = vn;
	vedge[edgei].endplane = planeB;
	edgeflag = 'N';
	break;
      }
      ++edgei;
    }
    if(edgeflag == 'Y') { // add edge
      add_vedge(vedge, edgenum, cont, planeA, planeC, planeB, poly, vn);
      ++edgenum;
    }

    // check edge (planeB, planeC)
    edgeflag = 'Y';
    edgei = startedge+1;
    while(edgei < edgenum) {
      if(((vedge[edgei].plane[0] == planeB)&&(vedge[edgei].plane[1] == planeC)) ||
	 ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeB))) {
	// already on list, add current vertex as endpt
	vedge[edgei].endpt = vn;
	vedge[edgei].endplane = planeA;
	edgeflag = 'N';
	break;
      }
      ++edgei;
    }
    if(edgeflag == 'Y') { // add edge
      add_vedge(vedge, edgenum, cont, planeB, planeC, planeA, poly, vn);
      ++edgenum;
      
      // ===== failsafe - if solution is not converging, perturb atom  =====
      // ===== coordinates and recalculate.                            =====
      if(edgenum >= 200) {
	//printf("********* invalid solution for hull, recalculating *********\n");
	seed[atomzero*3] = -1;  // reset to no seed vertex
	origcoor[0] = PDB[atomzero].coor[0];
	origcoor[1] = PDB[atomzero].coor[1];
	origcoor[2] = PDB[atomzero].coor[2];
	// perturb atom coordinates
	PDB[atomzero].coor[0] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
	PDB[atomzero].coor[1] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
	PDB[atomzero].coor[2] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
	recalc = 'Y';
	goto RESTART;
      }
    }
    ++vn;
  }

  /*--------------------------------------------------*/
  /*  now have voronoi polyhedron around given atom.  */
  /*  remove vertices outside of sphere, and          */
  /*  calculate intersection points with sphere.      */ 
  /*--------------------------------------------------*/

  // flag edges that may cross sphere boundary
  for(edgei=0; edgei<edgenum; ++edgei) {
    if((rado < poly[vedge[edgei].startpt].dist) || (rado < poly[vedge[edgei].endpt].dist)) {
      // one or both vertices fall outside of sphere
      arcflag = 'Y';
      vedge[edgei].arc = '?';
    } else {
      vedge[edgei].arc = 'X';
    }
  }
  
  // calculate new arc points
  for(edgei=0; edgei<edgenum; ++edgei) {
    if(vedge[edgei].arc != 'X') {
      if(solve_2xS(cont[vedge[edgei].plane[0]], cont[vedge[edgei].plane[1]], rado, arcpt0, arcpt1) == -1) {
	vedge[edgei].arc = 'X';  // mark edge as no associated arc point
	continue;
      }
      // test new arc points vs. adjacent planes, add if ok.
      testpA = vedge[edgei].startplane;
      testpB = vedge[edgei].endplane;
      testvalA = cont[testpA].Ai[0]*arcpt0[0] + cont[testpA].Ai[1]*arcpt0[1]
	+ cont[testpA].Ai[2]*arcpt0[2] + cont[testpA].Ai[3];
      testvalB = cont[testpB].Ai[0]*arcpt0[0] + cont[testpB].Ai[1]*arcpt0[1] 
	+ cont[testpB].Ai[2]*arcpt0[2] + cont[testpB].Ai[3];
      if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
	add_vertex(poly, vn, arcpt0, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
	poly[vn].dist = rado;
	++vn;
      }
      testvalA = cont[testpA].Ai[0]*arcpt1[0] + cont[testpA].Ai[1]*arcpt1[1] 
	+ cont[testpA].Ai[2]*arcpt1[2] + cont[testpA].Ai[3];
      testvalB = cont[testpB].Ai[0]*arcpt1[0] + cont[testpB].Ai[1]*arcpt1[1] 
	+ cont[testpB].Ai[2]*arcpt1[2] + cont[testpB].Ai[3];
      if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
	add_vertex(poly, vn, arcpt1, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
	poly[vn].dist = rado;
	++vn;
      }
    }
  }

  // reduce poly vertex list
  vj=0;
  for(vi=0; vi<vn; ++vi) {
    if(poly[vi].dist <= rado) {
      poly[vj] = poly[vi];
      ++vj;
    }
  }
  vn = vj;
 
  // ----- calculate center points and mark engulfing contacts -----
  for(cai=0; cai<NC; ++cai) {
    centerpt[cai].xi[0] = -cont[cai].Ai[0]*cont[cai].Ai[3];
    centerpt[cai].xi[1] = -cont[cai].Ai[1]*cont[cai].Ai[3];
    centerpt[cai].xi[2] = -cont[cai].Ai[2]*cont[cai].Ai[3];
    centerpt[cai].dist = fabs(cont[cai].Ai[3]);
    if(cont[cai].Ai[3] > 0.0) {
      cont[cai].flag = 'E';
    }
  }

  if(recalc == 'N') {
    save_seeds(cont, poly, vn, atomzero);
  } else {
    // reset atom coordinates to original values
    PDB[atomzero].coor[0] = origcoor[0];
    PDB[atomzero].coor[1] = origcoor[1];
    PDB[atomzero].coor[2] = origcoor[2];
  }
  return(vn);
}



/****************************
 * subroutine get_firstvert
 ****************************/

void get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC, int atomzero)
{
  long int seedi;
  int cai;
  double mindist;
  double ptA[3];
  double vectA[4];     // 4 values so it can be used as a plane equation as well
  double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
  double vtmin;        // minimum value for vt 
  double vtdiv;        // check for division by zero in vt calculation 
  double temppt[3];
  double ptdist;

  // get previous seed vertex if present
  seedi = atomzero*3;
  *planeA = -1;
  *planeB = -1;
  *planeC = -1;
  if(seed[seedi] != -1) {
    for(cai=0; cai<NC; ++cai) {
      if(cont[cai].index == seed[seedi]) {
	*planeA = cai;
      } else if(cont[cai].index == seed[seedi+1]) {
	*planeB = cai;
      } else if(cont[cai].index == seed[seedi+2]) {
	*planeC = cai;
      }
    }
  }
  if((*planeA != -1)&&(*planeB != -1)&&(*planeC != -1)) {
    return;
  } else {

    // -------------  find initial edge, on plane closest to origin  -------------
    mindist = 9.9e+9;  // dummy value
    for(cai=0; cai<NC; ++cai) {
      if(cont[cai].dist < mindist) {
	mindist = cont[cai].dist;
	*planeA = cai;
      }
    }

    mindist = 9.9e+9;  // dummy value
    for(cai=0; cai<NC+4; ++cai) {
      if(cai != *planeA) {
	vectA[0] = cont[*planeA].Ai[1]*cont[cai].Ai[2] - cont[*planeA].Ai[2]*cont[cai].Ai[1];
	vectA[1] = cont[*planeA].Ai[2]*cont[cai].Ai[0] - cont[*planeA].Ai[0]*cont[cai].Ai[2];
	vectA[2] = cont[*planeA].Ai[0]*cont[cai].Ai[1] - cont[*planeA].Ai[1]*cont[cai].Ai[0];
	vectA[3] = 0.0;
	if(solve_3x3(cont[*planeA].Ai, cont[cai].Ai, vectA, temppt) != -1) {
	  ptdist = sqrt(temppt[0]*temppt[0] + temppt[1]*temppt[1] + temppt[2]*temppt[2]);
	  if(ptdist < mindist) {
	    *planeB = cai;
	    mindist = ptdist;
	    ptA[0] = temppt[0];
	    ptA[1] = temppt[1];
	    ptA[2] = temppt[2];
 	  }
	}
      }
    }
    // recalc vector normal to planes A and B
    vectA[0] = cont[*planeA].Ai[1]*cont[*planeB].Ai[2] - cont[*planeA].Ai[2]*cont[*planeB].Ai[1];
    vectA[1] = cont[*planeA].Ai[2]*cont[*planeB].Ai[0] - cont[*planeA].Ai[0]*cont[*planeB].Ai[2];
    vectA[2] = cont[*planeA].Ai[0]*cont[*planeB].Ai[1] - cont[*planeA].Ai[1]*cont[*planeB].Ai[0];
    vectA[3] = 0.0;
  
    // get starting vertex on polyhedron
    vtmin = 9.9e+9;   // dummy value
    for(cai=0; cai<NC+4; ++cai) {
      if((cai != *planeA) && (cai != *planeB)) {
	vtdiv = (cont[cai].Ai[0]*vectA[0] +cont[cai].Ai[1]*vectA[1] +cont[cai].Ai[2]*vectA[2]);
	if(vtdiv != 0.0) {
	  vt = -(cont[cai].Ai[0]*ptA[0] +cont[cai].Ai[1]*ptA[1] +cont[cai].Ai[2]*ptA[2] +cont[cai].Ai[3])/vtdiv;
	  if(fabs(vt) < vtmin) {
	    vtmin = fabs(vt);
	    *planeC = cai;
	  }
	}
      }
    }
  }
  return;
}


/*************************
 * subroutine save_seeds
 *************************/

// saves starting vertices for atoms to be done

void save_seeds(struct plane cont[], struct vertex poly[], int NV, int atomzero)
{
  int vi;
  int seedi;

  for(vi=0; vi<NV; ++vi) {
    if(poly[vi].plane[2] != -1) {
      seedi = 3*cont[poly[vi].plane[0]].index;
      if(seed[seedi] == -1) {
	seed[seedi] = atomzero;
	seed[seedi+1] = cont[poly[vi].plane[1]].index;
	seed[seedi+2] = cont[poly[vi].plane[2]].index;
      }
      seedi = 3*cont[poly[vi].plane[1]].index;
      if(seed[seedi] == -1) {
	seed[seedi] = atomzero;
	seed[seedi+1] = cont[poly[vi].plane[0]].index;
	seed[seedi+2] = cont[poly[vi].plane[2]].index;
      }
      seedi = 3*cont[poly[vi].plane[2]].index;
      if(seed[seedi] == -1) {
	seed[seedi] = atomzero;
	seed[seedi+1] = cont[poly[vi].plane[0]].index;
	seed[seedi+2] = cont[poly[vi].plane[1]].index;
      }
    }
  }  
  return;
}



/*******************************
 * subroutine add_vedge
 *******************************/

// adds a new edge to the edge list. Direction of vector is tested so that
// it points away from the startpt, along the body of the polyhedron.

// stores information in variable vedge[edgenum].

void add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1, 
	       int testplane, struct vertex poly[], int startpt)
{
  double testpt[3];
  double testval;

  vedge[edgenum].V[0] = cont[plane0].Ai[1]*cont[plane1].Ai[2] - cont[plane0].Ai[2]*cont[plane1].Ai[1];
  vedge[edgenum].V[1] = cont[plane0].Ai[2]*cont[plane1].Ai[0] - cont[plane0].Ai[0]*cont[plane1].Ai[2];
  vedge[edgenum].V[2] = cont[plane0].Ai[0]*cont[plane1].Ai[1] - cont[plane0].Ai[1]*cont[plane1].Ai[0];
  vedge[edgenum].startpt = startpt;
  vedge[edgenum].endpt = -1;  // flag, edge not completed.
  vedge[edgenum].plane[0] = plane0;
  vedge[edgenum].plane[1] = plane1;
  vedge[edgenum].startplane = testplane;
  vedge[edgenum].arc = '.'; // dummy value.

  // test direction of vector
  testpt[0] = poly[startpt].xi[0] + vedge[edgenum].V[0];
  testpt[1] = poly[startpt].xi[1] + vedge[edgenum].V[1];
  testpt[2] = poly[startpt].xi[2] + vedge[edgenum].V[2];
  
  testval = cont[testplane].Ai[0]*testpt[0] +cont[testplane].Ai[1]*testpt[1] 
    + cont[testplane].Ai[2]*testpt[2] +cont[testplane].Ai[3];
  if(testval > 0.0) { // vector is in wrong direction
    vedge[edgenum].V[0] = -vedge[edgenum].V[0];
    vedge[edgenum].V[1] = -vedge[edgenum].V[1];
    vedge[edgenum].V[2] = -vedge[edgenum].V[2];
  }
  return;
}

/*******************************
 * subroutine add_vertex
 *******************************/

// add polyhedron vertex to local list for current atom contacts

int add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC)
{ 
  poly[vn].xi[0] = coor[0];
  poly[vn].xi[1] = coor[1];
  poly[vn].xi[2] = coor[2];
  poly[vn].plane[0] = planeA;
  poly[vn].plane[1] = planeB;
  poly[vn].plane[2] = planeC;
  poly[vn].dist = sqrt(coor[0]*coor[0] +coor[1]*coor[1] +coor[2]*coor[2]);
  return(0);
}

/*********************************
 * function order_faces          
 *********************************/

// output is array ptorder, the order of points around each face.
// return values are 'I' for internal atom, 'S' for surface atom.

char order_faces(int atomzero,struct vertex poly[], struct vertex centerpt[], float rado,
		 int NC, int NV, struct plane cont[], struct ptindex ptorder[])
{
  int planeX;             // current plane, 1 to NC
  int planeY;             // second plane, to get adjacent points on polygon
  int surfcount;          // number of points defining a given plane of contact
  int vi, vi2;            // vertices counter
  int tempsi;             // temporary storage for exchanging indices
  double tempcos;         // for exchanging cosines
  double cos10X[50];      // cos of angle between points poly[1],poly[0],and poly[X]
  double temppt[3];       // temp coordinates of new arc point
  char surfatom;          // return value: 'I' internal atom, 'S' surface atom

  surfatom = 'I'; // internal atom
  for(vi=0; vi<NV; ++vi) {
    if(poly[vi].plane[2] == -1) {
      surfatom = 'S'; // surface atom 
      break;
    }
  }
  // for surface calculation only
  // if(surfatom == 'I') return('I');
  
  for(planeX=0; planeX < NC; ++planeX) {
    if(cont[planeX].flag == 'X') { // hidden
      ptorder[planeX].numpts = 0;
      continue;
    }

    surfcount=0;
    for(vi=0; vi<NV; ++vi) {
      // index all points comprising surface for planeX
      if((poly[vi].plane[0]==planeX) || (poly[vi].plane[1]==planeX) || (poly[vi].plane[2]==planeX)) {
	ptorder[planeX].pt[surfcount] = vi;
	++surfcount;
      }
    }

    if(surfcount > 2) {

      /*-------------------------------------------------------------------*/
      /* three or more surface points, need additional point so no arcs    */
      /* are greater than 180 deg. Make point opposite first pt on an arc. */
      /* (if no points on arc, extra point isn't needed.)                  */
      /*-------------------------------------------------------------------*/

      for(vi=0; vi<surfcount; ++vi) {
	if(poly[ptorder[planeX].pt[vi]].plane[2] == -1) { // on arc, calc pt. opposite

	  // get coordinates of point
	  temppt[0] = 2.0*centerpt[planeX].xi[0] - poly[ptorder[planeX].pt[vi]].xi[0];
	  temppt[1] = 2.0*centerpt[planeX].xi[1] - poly[ptorder[planeX].pt[vi]].xi[1];
	  temppt[2] = 2.0*centerpt[planeX].xi[2] - poly[ptorder[planeX].pt[vi]].xi[2];
	  
	  //keep point if it's valid
	  if(test_point(temppt, cont, NC, rado, planeX, -1, -1) == 'Y') {
	    poly[NV].xi[0] = temppt[0];
	    poly[NV].xi[1] = temppt[1];
	    poly[NV].xi[2] = temppt[2];
	    poly[NV].plane[0] = planeX;
	    poly[NV].plane[1] = -1;
	    poly[NV].plane[2] = -1;
	    poly[NV].dist = rado;
	    ptorder[planeX].pt[surfcount] = NV;
	    ++surfcount;
	    ++NV;
	  }
	  break;
	}
      }
    }
    ptorder[planeX].numpts = surfcount;

    if(surfcount > 3) {
      // get two points on same line (two common planes).
      // all points already share one plane (planeX), find another.
      if(poly[ptorder[planeX].pt[0]].plane[0] == planeX) {
	planeY = poly[ptorder[planeX].pt[0]].plane[1];
      } else {
	planeY = poly[ptorder[planeX].pt[0]].plane[0];
      }
      
      //find another point on the same line (2 common planes)
      for(vi=1; vi<surfcount; ++vi) {
	if((poly[ptorder[planeX].pt[vi]].plane[0]==planeY) || (poly[ptorder[planeX].pt[vi]].plane[1]==planeY) 
	   || (poly[ptorder[planeX].pt[vi]].plane[2]==planeY)) {
	  break;
	}
      }    

      //swap index for pt[1] and pt[vi], so points 0 and 1 are on same line
      tempsi = ptorder[planeX].pt[vi];
      ptorder[planeX].pt[vi] = ptorder[planeX].pt[1];
      ptorder[planeX].pt[1] = tempsi;

      // calculate cosine between points indexed 1,0,X
      for(vi=2; vi<surfcount; ++vi) {
	cos10X[vi] = cosPQR(poly[ptorder[planeX].pt[1]].xi, poly[ptorder[planeX].pt[0]].xi, 
			    poly[ptorder[planeX].pt[vi]].xi);
      }

      // order by cosines, decreasing order
      for(vi=2; vi<surfcount-1; ++vi) {
	for(vi2=vi+1; vi2<surfcount; ++vi2) {
	  if(cos10X[vi] < cos10X[vi2]) {
	    // swap indices if points in wrong order
	    tempsi = ptorder[planeX].pt[vi];
	    ptorder[planeX].pt[vi] = ptorder[planeX].pt[vi2];
	    ptorder[planeX].pt[vi2] = tempsi;
	    tempcos = cos10X[vi];
	    cos10X[vi] = cos10X[vi2];
	    cos10X[vi2] = tempcos;
	  }
	}
      }
    }
  }
  return(surfatom);
}


/**************************
 * subroutine calc_areas
 **************************/

void calc_areas(struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, 
		struct plane cont[], struct ptindex ptorder[], int atomzero)
{
  char   engflag;         // ='Y' if an engulfing plane is present
  int    planeX;          // current plane, 1 to NC
  int    NP;              // number of points on face
  int    vi;              // vertices counter
  double area;            // area of polygon
  int    pa1, pa2;        // possible common planes
  int    commplane;       // plane shared by adjacent points (not planeX)
  int    epi;             // engulfing plane counter
  int    engplane[4];     // index to engulfed planes
  double maxSAS;          // maximum solvent exposed surface, no atoms other than engulfing.
  struct vertex ptB, ptC; // arc point intersections
  int    currpt, nextpt;
  double cosNN1[40];      // angle between vertex N and vertex N+1 
  double cosNzero[40];    // angle between vertex N and vertex N+1 
  double tanprod;         // product of tangents
  int    v0, va, vb;      // vertices for arc calculation

  double U,V,W,X;
  double tansqrS, tansqrSA, tansqrSB, tansqrSC;

  engflag = 'N';
  epi = 0;

  //RESET AREAS TO ZERO
  for(planeX=0; planeX < NC; ++planeX) {
    cont[planeX].area = 0.0;
    if(cont[planeX].flag == 'E') {
      engflag = 'Y';
      engplane[epi] = planeX;
      ++epi;
    }
  }

  if(engflag == 'Y') { // engulfing plane correction - project points onto sphere surface.
    project_points(poly, centerpt, rado, NC, NV, cont);
  }

  /* ---------------------------- */
  /* calculate area for each face */
  /* ---------------------------- */

  for(planeX=0; planeX < NC; ++planeX) {
    NP = ptorder[planeX].numpts;
    area = 0.0;
    if(cont[planeX].flag == 'X') {
      continue;
    }

    // if there are no points on a valid contact, area is spherical cap
    if(NP == 0) {
      if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') {
	cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist);
      }
    } else if(NP == 2) {  // only two contact points, check which part of arc
      if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') { // area is (cap - arc)
	cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist)
	  - spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]], 
			  poly[ptorder[planeX].pt[1]], rado);
      } else {            // area is arc.
	cont[planeX].area = spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]], 
					  poly[ptorder[planeX].pt[1]], rado);
      }
    } else {

      // ------ calculate cosines and angles ------
      for(vi=0; vi<NP; ++vi) {
	v0 = ptorder[planeX].pt[0];
	va = ptorder[planeX].pt[vi];
	vb = ptorder[planeX].pt[(vi+1)%NP];

	// calculate cosines between adjacent vertices
	cosNN1[vi] = (( poly[va].xi[0]*poly[vb].xi[0] + poly[va].xi[1]*poly[vb].xi[1] 
			+ poly[va].xi[2]*poly[vb].xi[2] ) / (poly[va].dist*poly[vb].dist));

	// calculate cosines between vertex zero and vertex 'vi'
	if(vi != 0) {
	  cosNzero[vi] = ((poly[v0].xi[0]*poly[va].xi[0] + poly[v0].xi[1]*poly[va].xi[1] 
			   + poly[v0].xi[2]*poly[va].xi[2] ) / (poly[v0].dist*poly[va].dist));
	}
      }
      
      // ----- calculate area of triangles in face -----
      for(vi=1; vi<(NP-1); ++vi) {
	U = sqrt((1+cosNzero[vi])*(1+cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
	V = sqrt((1-cosNzero[vi])*(1-cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
	W = sqrt((1-cosNzero[vi])*(1+cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
	X = sqrt((1+cosNzero[vi])*(1-cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
	tansqrS  = (1-U+V+W+X)/(1+U-V-W-X);
	tansqrSA = (1-U-V-W+X)/(1+U+V+W-X);
	tansqrSB = (1-U-V+W-X)/(1+U+V-W+X);
	tansqrSC = (1-U+V-W-X)/(1+U-V+W+X);
	tanprod = sqrt(tansqrS*tansqrSA*tansqrSB*tansqrSC);
	if(tanprod > 0.0) {
	  area += 4.0*rado*rado*atan(sqrt(tanprod));
	}
      }

      // ----- add area of arc segments  -----
      for(vi=0; vi<NP; ++vi) {
	va = ptorder[planeX].pt[vi];
	vb = ptorder[planeX].pt[(vi+1)%NP];
	
	//check if adjacent points are arc segments
	if((poly[va].plane[2] == -1) && (poly[vb].plane[2] == -1)) {
	  // if on same two planes, no arc.
	  if((poly[va].plane[0]+poly[va].plane[1]) != (poly[vb].plane[0]+poly[vb].plane[1])) {
	    area += spherical_arc(centerpt[planeX], poly[va], poly[vb], rado);
	  }
	}
      }
      cont[planeX].area = area;
    }
  }

  // --------------------------------------------------------
  //  add correction terms for engulfing planes, if required
  // --------------------------------------------------------

  if(engflag == 'Y') {
    for(planeX=0; planeX < NC; ++planeX) {
      if(cont[planeX].flag != 'E') {
	continue;
      }

      NP = ptorder[planeX].numpts;
      for(vi=0; vi<NP; ++vi) {
	currpt = ptorder[planeX].pt[vi];
	nextpt = ptorder[planeX].pt[(vi+1)%NP];

	// find common second plane, if any.
	if(poly[currpt].plane[0] == planeX) {
	  pa1 = poly[currpt].plane[1];
	  pa2 = poly[currpt].plane[2];
	} else {
	  pa1 = poly[currpt].plane[0];
	  if(poly[currpt].plane[1] == planeX) {
	    pa2 = poly[currpt].plane[2];
	  } else {
	    pa2 = poly[currpt].plane[1];
	  }
	}

	if((pa1 == poly[nextpt].plane[0]) || (pa1 == poly[nextpt].plane[1])
	   || (pa1 == poly[nextpt].plane[2])) {
	  commplane = pa1;
	} else if((pa2 == poly[nextpt].plane[0]) || (pa2 == poly[nextpt].plane[1])
		  || (pa2 == poly[nextpt].plane[2])) { 
	  commplane = pa2;
	} else {
	  continue;
	}
	if((commplane != -1) && (cont[commplane].flag != 'E')) {
	  // add correction to commplane area. here centerpt is from engulfing plane.
	  cont[commplane].area += spherical_arc(centerpt[planeX], poly[currpt], poly[nextpt], rado);
	  if(NP == 2) break;  // otherwise would repeat adding area
	}
      }
    }

    // -----------------------------------------------
    // ------ calculate engulfed contact areas -------
    // -----------------------------------------------

    if(epi == 1) {
      cont[engplane[0]].area = 2.0*PI*rado*(rado+cont[engplane[0]].dist);
    } else if(epi == 2) { 
      if(solve_2xS(cont[engplane[0]],cont[engplane[1]], rado, ptB.xi, ptC.xi)== -1) {
	cont[engplane[0]].area = 2.0*PI*rado*rado;
	cont[engplane[1]].area = 2.0*PI*rado*rado;
      } else {
	ptB.dist = rado;
	ptC.dist = rado;
	maxSAS = spherical_arc(centerpt[engplane[0]], ptB, ptC, rado);
 	maxSAS += spherical_arc(centerpt[engplane[1]], ptB, ptC, rado);
	cont[engplane[0]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
	cont[engplane[1]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
      }
    } else if(epi>=3) {
      // no exposed surface if there are three or more engulfing contacts 
      for(planeX=0; planeX<NC; ++planeX) {
	if(cont[planeX].flag == 'E') {
	  cont[planeX].area = 4.0*PI*rado*rado/epi;
	} else {
	  cont[planeX].area = 0.0;
	}
      }
    }
  }
  return;
}


/***************************
 * subroutine save_areas
 ***************************/

void save_areas(struct plane cont[], struct contactlist contlist[], int NC, int atomzero)
{
  int cai;
  long int currindex, previndex, nextindex;

  if(numcarec > (ca_recsize-100)) {
    ca_recsize += 10000;
    ca_rec = realloc(ca_rec, ca_recsize*sizeof(struct ca_struct));
    if(!ca_rec) { 
      printf("memory allocation error\n"); 
      exit(1);
    }
  }

  // first overwrite previous records
  cai = 0;
  currindex = ca_index[atomzero];
  previndex = -1;
  while((currindex != -1) && (cai<NC)) {
    if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
      ca_rec[currindex].atom = cont[cai].index;
      ca_rec[currindex].area = cont[cai].area;
      ca_rec[currindex].dist = contlist[cai].dist;
      nextindex = ca_rec[currindex].prev; // next index is prev record from old list
      ca_rec[currindex].prev = previndex;
      previndex = currindex;
      currindex = nextindex;
    }
    ++cai;
  }
  ca_index[atomzero] = previndex;
  
  // then add new records
  while(cai<NC) {
    if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
      ca_rec[numcarec].atom = cont[cai].index;
      ca_rec[numcarec].area = cont[cai].area;
      ca_rec[numcarec].dist = contlist[cai].dist;
      ca_rec[numcarec].prev = ca_index[atomzero];
      ca_index[atomzero] = numcarec;
      ++numcarec;
    }
    ++cai;
  }
  PDB[atomzero].done = 'Y';

  // index contacts for atoms not yet done
  for(cai=0; cai<NC; ++cai) {
    if(PDB[cont[cai].index].done != 'Y') {
      if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
	ca_rec[numcarec].atom = atomzero;
	ca_rec[numcarec].prev = ca_index[cont[cai].index];
	ca_index[cont[cai].index] = numcarec;
	++numcarec;
      }
    }
  }
  return;
}

/*****************************
 * subroutine project_points
 *****************************/

// this subroutine corrects for engulfing atoms. The center of the engulfing plane
// contact is used as the center of projection instead of the center of the atom.
// points already on the surface are not moved, preserving the SAS.

void project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC, 
		    int NV, struct plane cont[])
{
  int epi;               // engulfing plane counter
  int cai;               // contact atom counter
  int engplane[4];       // index to engulfing planes
  double projpt[3];      // center point for projecting internal pts onto surface of sphere
  double pt0[3], pt1[3]; // temporary points for intersection solutions
  double V[3];           // vector from projection point to vertex
  double *P;             // pointer for vertex coordinates to be projected
  double a,b,c,k;        // for solving quadratic eqn for point projection
  int vi;                // vertex counter

  // count and mark engulfing planes
  epi = 0;
  for(cai=0; cai<NC; ++cai) {
    if(cont[cai].flag == 'E') {
      engplane[epi] = cai;
      ++epi;
    }
  }

  // get projpt[] for projecting points to surface
  if(epi == 1) {
    projpt[0] = centerpt[engplane[0]].xi[0];
    projpt[1] = centerpt[engplane[0]].xi[1];
    projpt[2] = centerpt[engplane[0]].xi[2];
  } else if(epi == 2) {
    solve_2xS(cont[engplane[0]],cont[engplane[1]],rado, pt0, pt1);
    projpt[0] = (pt0[0]+pt1[0])/2;
    projpt[1] = (pt0[1]+pt1[1])/2;
    projpt[2] = (pt0[2]+pt1[2])/2;
  } else {
    solve_3x3(cont[engplane[0]].Ai, cont[engplane[1]].Ai, cont[engplane[2]].Ai, pt0);
    projpt[0] = pt0[0];
    projpt[1] = pt0[1];
    projpt[2] = pt0[2];
  }

  for(vi=0; vi<NV; ++vi) {
    if(poly[vi].plane[2] != -1) {
      // project point to intersection of engulfing plane(s) and surface of sphere
      P = poly[vi].xi;
      V[0] = P[0] - projpt[0];
      V[1] = P[1] - projpt[1];
      V[2] = P[2] - projpt[2];
      a = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
      b = 2*(P[0]*V[0] +P[1]*V[1] +P[2]*V[2]);
      c = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] - rado*rado;  // c is < 0 
      k = (sqrt(b*b - 4.0*a*c) - b)/(2*a);                // k is > 0
      P[0] += k*V[0];
      P[1] += k*V[1];
      P[2] += k*V[2];      
      poly[vi].dist = rado;
    }
  }
  return;
}


/************************
 * function test_point     
 ************************/

// this function tests a given point versus a set of plane constraints
// it returns a value of 'Y' if the point is OK, and 'N' if there
// was a violation of the plane inequality.
// ptX =(xo, yo, zo); if Axo+Byo+Czo+D > 0, point is behind plane (hidden)
// planes A,B,C are the planes that the point lies on, don't test.

char test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC)
{
  int cp;      // counter, number of planes

  // if point is not behind any plane, keep point.
  for(cp=0; cp<NC; ++cp) {
    //if pt is on current plane, get next plane
    if((cp != planeA) && (cp != planeB) && (cp != planeC) && (cont[cp].flag != 'X')) {
      if((cont[cp].Ai[0]*ptX[0] + cont[cp].Ai[1]*ptX[1] + cont[cp].Ai[2]*ptX[2] + cont[cp].Ai[3]) > 0.0) {
	//point is behind plane, not on polyhedron.
	return('N');
      }
    }
  }
  return('Y');
}


/****************************
 * function solve_3x3   
 ****************************/

// determines the intersection of three planes
// (solves a system of three linear equations and 3 unknowns)
// input is 3 four element arrays, representing Ax+By+Cz+D=0
// output is a three element array, (xi,xj,xk).

int solve_3x3(double eq0[], double eq1[], double eq2[], double pt[])
{
  double cof00, cof01, cof02;   // matrix cofactors
  double cof10, cof11, cof12;
  double cof20, cof21, cof22;
  double det;                   // determinant of matrix

  cof00 =  eq1[1]*eq2[2] - eq2[1]*eq1[2];
  cof01 = -eq1[0]*eq2[2] + eq2[0]*eq1[2];
  cof02 =  eq1[0]*eq2[1] - eq2[0]*eq1[1];

  cof10 = -eq0[1]*eq2[2] + eq2[1]*eq0[2];
  cof11 =  eq0[0]*eq2[2] - eq2[0]*eq0[2];
  cof12 = -eq0[0]*eq2[1] + eq2[0]*eq0[1];

  cof20 =  eq0[1]*eq1[2] - eq1[1]*eq0[2];
  cof21 = -eq0[0]*eq1[2] + eq1[0]*eq0[2];
  cof22 =  eq0[0]*eq1[1] - eq1[0]*eq0[1];

  det = eq0[0]*cof00 + eq0[1]*cof01 + eq0[2]*cof02;
  if(det == 0.0) {
    //printf("no solution for equation set, determinant is zero\n");
    return(-1);
  } else {
    pt[0] = -(eq0[3]*cof00 + eq1[3]*cof10 + eq2[3]*cof20)/det;
    pt[1] = -(eq0[3]*cof01 + eq1[3]*cof11 + eq2[3]*cof21)/det;
    pt[2] = -(eq0[3]*cof02 + eq1[3]*cof12 + eq2[3]*cof22)/det;
    return(0);
  }
}

/**********************************
 * Function solve_2xS       
 * revised 19/02/2001  BJM
 ***********************************/

// determines the intersection of two planes and a sphere radius 'rad'
// input is 2 four element arrays, representing Ax+By+Cz+D=0
// plus the radius of a sphere centered on the origin.
// output is two three element arrays pt0 and pt1, (xi,xj,xk).
// return value is -1 if no real solution exists.

int solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[])
{
  double eq2[3];              // eqn of plane through (0,0,0) and perp. to other two
  double cof00, cof01, cof02; // matrix cofactors
  double cof10, cof11, cof12; // (don't need cof20, cof21, cof22)
  double det;                 // determinant of matrix
  double avgpt[3];            // average of two solution points
  double t;                   // parameter in eqn of line: x=xo+At, y=yo+Bt, z=zo+Ct. 
  int xi;                     // coordinate counter

  eq2[0] = eq0.Ai[1]*eq1.Ai[2] - eq0.Ai[2]*eq1.Ai[1];
  eq2[1] = eq0.Ai[2]*eq1.Ai[0] - eq0.Ai[0]*eq1.Ai[2];
  eq2[2] = eq0.Ai[0]*eq1.Ai[1] - eq0.Ai[1]*eq1.Ai[0];

  cof00 =  eq1.Ai[1]*eq2[2] - eq2[1]*eq1.Ai[2];
  cof01 = -eq1.Ai[0]*eq2[2] + eq2[0]*eq1.Ai[2];
  cof02 =  eq1.Ai[0]*eq2[1] - eq2[0]*eq1.Ai[1];

  cof10 = -eq0.Ai[1]*eq2[2] + eq2[1]*eq0.Ai[2];
  cof11 =  eq0.Ai[0]*eq2[2] - eq2[0]*eq0.Ai[2];
  cof12 = -eq0.Ai[0]*eq2[1] + eq2[0]*eq0.Ai[1];

  det = eq2[0]*eq2[0] + eq2[1]*eq2[1] + eq2[2]*eq2[2];
  if(det == 0) {
    //printf("no solution in solve_2xS\n");
    return(-1);
  }

  avgpt[0] = -(eq0.Ai[3]*cof00 + eq1.Ai[3]*cof10)/det;
  avgpt[1] = -(eq0.Ai[3]*cof01 + eq1.Ai[3]*cof11)/det;
  avgpt[2] = -(eq0.Ai[3]*cof02 + eq1.Ai[3]*cof12)/det;

  t = (rado*rado-avgpt[0]*avgpt[0]-avgpt[1]*avgpt[1]-avgpt[2]*avgpt[2])/det;
  if(t<0) {
    return(-1);
  } else {
    t = sqrt(t);
  }
  for(xi=0; xi<3; ++xi) {
    pt0[xi] = avgpt[xi] + t*eq2[xi];
    pt1[xi] = avgpt[xi] - t*eq2[xi];
  }
  return(0);
}


/***************************
 * function cosPQR         *
 ***************************/

// this function returns the cosine of the angle between three points P,Q,R 
// with Q being the center point.

double cosPQR( double ptP[], double ptQ[], double ptR[])
{
  double QP[3];    // vector from Q to P
  double QR[3];    // vector from Q to R
  double cosine;   // cosine of angle PQR at Q.

  // calculate vectors
  QP[0] = ptP[0] - ptQ[0];
  QP[1] = ptP[1] - ptQ[1];
  QP[2] = ptP[2] - ptQ[2];
  QR[0] = ptR[0] - ptQ[0]; 
  QR[1] = ptR[1] - ptQ[1]; 
  QR[2] = ptR[2] - ptQ[2]; 

  //calculate cosine
  cosine = (QP[0]*QR[0]+ QP[1]*QR[1]+ QP[2]*QR[2])
    /sqrt((QP[0]*QP[0]+ QP[1]*QP[1]+ QP[2]*QP[2]) * (QR[0]*QR[0]+ QR[1]*QR[1]+ QR[2]*QR[2]));

  return(cosine);
}


/********************************
 * function spherical_arc       *	   
 ********************************/

// given two points in cartesian space and the center of the spherical
// cap between atom A and the origin, plus the radius of the sphere,
// this function returns the area of an arc between points B and C.
// the sides of the arc are the great circle distance (shortest distance)
// and the arc of the spherical cap centered on line AO.

double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado)
{
  // here, cosAOC=cosAOB, AOB = AOC.
  // BAC is the angle at the point on the origin, on line AO
  
  double cosAOB, cosBOC; // cosines of angles centered on (0,0,0)
  double cosBAC, angBAC;    // angle and cosine at vertex of cap
  double U, V, W;
  double tansqrS, tansqrSA, tansqrSB;
  double tanprod;        // product of tangents for spherical triangle
  double area;           // the value to be returned

  cosAOB = (ptAo.xi[0]*ptB.xi[0] + ptAo.xi[1]*ptB.xi[1] + ptAo.xi[2]*ptB.xi[2])/(ptAo.dist*ptB.dist);
  cosBOC = (ptB.xi[0]*ptC.xi[0] + ptB.xi[1]*ptC.xi[1] + ptB.xi[2]*ptC.xi[2])/(ptB.dist*ptC.dist);

  U = (1+cosAOB)*sqrt((1+cosBOC)/8.0);
  V = (1-cosAOB)*sqrt((1+cosBOC)/8.0);
  W = sqrt((1-cosAOB)*(1+cosAOB)*(1-cosBOC)/8.0);  // W == X
  tansqrS  = (1-U+V+W+W)/(1+U-V-W-W);
  tansqrSB = (1-U-V)/(1+U+V);
  tansqrSA = (1-U+V-W-W)/(1+U-V+W+W);
  if((tansqrS*tansqrSA) > 0.0) {
    tanprod = sqrt(tansqrSB*tansqrSB*tansqrS*tansqrSA);
  } else {
    tanprod = 0.0;
  }

  cosBAC = cosPQR(ptB.xi, ptAo.xi, ptC.xi);
  if(cosBAC>1.0) { 
    angBAC = 0.0; 
  } else if(cosBAC<-1.0) {
    angBAC = PI; 
  } else { 
    angBAC = acos(cosBAC);
  } 

  // area is area of spherical cap segment minus area of spherical triangle
  area = rado*(angBAC*(rado-ptAo.dist) - 4.0*rado*atan(sqrt(tanprod)));
  return(area);
}


/******************************
 *  function read_PDB         *
 ******************************/

// hydrogen atom records are skipped
// PDB[] is now a globally defined pointer.
// function fills array 'PDB' with data read from PDB file.
// return value is a pointer to the location of the PDB array.
// reads only one of multiple possible structures- if ATOM AltLoc
// record (column 17) is different than ' ' or 'A', atom is omitted. 

struct atom *read_PDB(char *FILE_ptr, struct atom *PDB, long int *numPDBatoms)
{
  FILE  *infile;            // input file pointer
  char  *oneline_ptr;       // to test for end of file
  char   oneline[89];       // the data read
  long int PDBtot;          // number of atoms read from the array
  long int PDBarraysize;

  // ------- open file for input ---------
  infile = fopen(FILE_ptr, "r");
  if (infile == NULL) {
    printf("Can't open file %s for reading\n", FILE_ptr);
    exit(8);
  }
  
  PDBarraysize = 1000;  // initial memory allocation
  PDB = malloc(PDBarraysize*sizeof(struct atom));
  if(!PDB) {
    printf("memory allocation for PDB atoms failed\n");
    exit(1);
  }

  PDBtot = 0;
  while (1) {
    oneline_ptr = fgets(oneline, sizeof(oneline), infile);
    if ((strncmp(oneline, "END",3)==0) || (oneline_ptr == NULL)) {
      fclose(infile);
      *numPDBatoms = PDBtot;
      assign_radii(PDB, PDBtot);
      return(PDB);
    }

    // if PDB array is getting full, add more elements to array
    if(PDBtot>=PDBarraysize) {
      PDBarraysize += 1000;
      PDB = realloc(PDB, PDBarraysize*sizeof(struct atom));
      if(!PDB) {
	printf("memory allocation for PDB atoms failed\n");
	exit(1);
      }
    }

    //**********include this line to include HETATMs -->
    if ((!strncmp(&oneline[0],"ATOM  ",6))||(!strncmp(&oneline[0], "HETATM",6))) {
    //**********include this line if only ATOMs -->
    //if (!strncmp(&oneline[0], "ATOM  ", 6)) {

      // ------- skip alternate location records ----------
      if ((oneline[16] != ' ') && (oneline[16] != 'A')) continue;
      // ------- skip HOH records ----------
      if (!strncmp(&oneline[17], "HOH", 3)) continue;
      // ------- skip explicit hydrogens -------
      if(oneline[13] == 'H') continue;

      PDB[PDBtot].atomnum = atoi(&oneline[6]);        // atom number
      strncpy(PDB[PDBtot].atomname, &oneline[12], 4); // atom name
      PDB[PDBtot].atomname[4] = '\0';
      strncpy(PDB[PDBtot].res, &oneline[17], 3);      // residue name
      PDB[PDBtot].res[3] = '\0';
      PDB[PDBtot].chain = oneline[21];
      PDB[PDBtot].resnum  = atoi(&oneline[22]);       // residue number
      PDB[PDBtot].coor[0] = atof(&oneline[30]);       // x coordinate
      PDB[PDBtot].coor[1] = atof(&oneline[38]);       // y coordinate
      PDB[PDBtot].coor[2] = atof(&oneline[46]);       // z coordinate
      ++PDBtot;
    }
  }   
  return(PDB);   // dummy return
}

/********************************
 * function assign_radii
 ********************************/

// uses radii from file radii.dat.
// edited to assign radii for RNA and DNA, 25-Sept-2003  BJM
// radius for Phosphate is same as default radii, using default
	 
void assign_radii(struct atom *PDB, long int PDBtot)
{
	long int atomi;
	// FILE  *radfile;            // radius data file name
	// char   filename[20] = "radii.dat";
	// char   typename[11][5] = {"C3H0","C3H1","C4  ","N3H0","N3H1","N3H2","N4  ","O1H0","O2H1","S   ","DEF "};
	/* ---------------------   0      1      2      3      4      5      6      7      8      9      10   */
	float  radius[11];
	// int    linenum= 0;
	// int    recnum = 0;
	// char  *oneline_ptr;       // to test for end of file
	// char   oneline[100];      // the data read

	// ------- open file for input ---------
	// radfile = fopen(filename, "r");
	// if (radfile == NULL) {
	// 	printf("Can't open file %s for reading\n", filename);
	// 	exit(8);
	// }
	// while(1) {
	// 	oneline_ptr = fgets(oneline, sizeof(oneline), radfile);
	// 	if (!oneline_ptr) {
	// 		fclose(radfile);
	// 		break;
	// 	}
	// 	if(oneline[0] != '#') {
	// 		if(strncmp(&oneline[2], typename[recnum], 4)) {
	// 			printf("in file '%s', expected atom type %s, line %d\n", filename, typename[recnum], linenum);
	// 			exit(1);
	// 		} else {
	// 			radius[recnum] = atof(&oneline[8]);
	// 			++recnum;
	// 		}
	// 	}
	// 	++linenum;
	// }
  radius[0] = 1.61; // C3H0
  radius[1] = 1.76; // C3H1
  radius[2] = 1.88; // C4
  radius[3] = 1.64; // N3H0
  radius[4] = 1.64; // N3H1
  radius[5] = 1.64; // N3H2
  radius[6] = 1.64; // N4
  radius[7] = 1.42; // O1H0
  radius[8] = 1.46; // O2H1
  radius[9] = 1.77; // S
  radius[10] = 1.80; // DEF



	for(atomi=0; atomi<PDBtot; ++atomi) {
		// ===============================================================
		// check if residue is DNA or RNA - ADDED  25-Sept-2003 BJM
		if(!strncmp(PDB[atomi].res, "  C", 3) || !strncmp(PDB[atomi].res, "  G", 3)
	 || !strncmp(PDB[atomi].res, "  A", 3) || !strncmp(PDB[atomi].res, "  T", 3)
	 || !strncmp(PDB[atomi].res, "  U", 3)) {
			// residue is DNA or RNA
			if(PDB[atomi].atomname[1] == 'P') {
				// PHOSPHORUS atom
				PDB[atomi].radius = radius[10];  // radius is same as default
			} else if(PDB[atomi].atomname[1] == 'O') {
				// OXYGEN atom
				// O2H0 is same radius as O2H1
				if(PDB[atomi].atomname[3] == 'P') { // part of phosphate
					PDB[atomi].radius = radius[7];  // O1H0
				} else if(PDB[atomi].atomname[3] == '*') { // ribose
					PDB[atomi].radius = radius[8];  // O2H1 ( also O2H0, same radius)
				} else if(PDB[atomi].atomname[3] == ' ') { // base
					PDB[atomi].radius = radius[7];  // O1H0
				}
			} else if(PDB[atomi].atomname[1] == 'N') {
				// NITROGEN ATOM
				if((PDB[atomi].atomname[2] == '2') || (PDB[atomi].atomname[2] == '4') || (PDB[atomi].atomname[2] == '6')) {
					PDB[atomi].radius = radius[5];  // N3H2
				} else { // part of ring
					PDB[atomi].radius = radius[3];  // N3H0 (== N2H0, same radius)
				}
			} else if(PDB[atomi].atomname[1] == 'C') {
				// CARBON ATOM
				if(PDB[atomi].atomname[3] == '*') { // ribose carbon
					PDB[atomi].radius = radius[2];  // C4
				} else if(PDB[atomi].atomname[2] == '2') { 
					if(PDB[atomi].res[2] == 'A') {
						PDB[atomi].radius = radius[1];  // C3H1
					} else {
						PDB[atomi].radius = radius[0];  // C3H0						
					}
				} else if(PDB[atomi].atomname[2] == '4') { 
					PDB[atomi].radius = radius[0];  // C3H0
				} else if(PDB[atomi].atomname[2] == '5') { 
					if((PDB[atomi].res[2] == 'C') || (PDB[atomi].res[2] == 'U')) {
						PDB[atomi].radius = radius[1];  // C3H1
					} else {
						PDB[atomi].radius = radius[0];  // C3H0
					}
				} else if(PDB[atomi].atomname[2] == '6') { 
					if((PDB[atomi].res[2] == 'A') || (PDB[atomi].res[2] == 'G')) {
						PDB[atomi].radius = radius[0];  // C3H0
					} else {
						PDB[atomi].radius = radius[1];  // C3H1
					}
				} else if(PDB[atomi].atomname[2] == '8') { 
					PDB[atomi].radius = radius[1];  // C3H1
				} else { 
					PDB[atomi].radius = radius[1];  // C3H1, dummy value					
				}
			}
			// ============================ end addition =========================
		} else {
			// if not DNA/RNA, then amino acid residue
			if(PDB[atomi].atomname[1] == 'O') {
				//   OXYGEN
				if((PDB[atomi].atomname[2] == 'G')||(PDB[atomi].atomname[2] == 'H')) {
					PDB[atomi].radius = radius[8];  // O2H1
				} else {
					PDB[atomi].radius = radius[7];  // O2H0
				}
			} else if(PDB[atomi].atomname[1] == 'S') {
				// SULFUR
				PDB[atomi].radius = radius[9]; // S
			} else if(PDB[atomi].atomname[1] == 'N') {
				// NITROGEN
				if(!strncmp(PDB[atomi].res, "PRO", 3)) {
					PDB[atomi].radius = radius[3];  // N3H0
				} else if(PDB[atomi].atomname[2] == ' ') {
					PDB[atomi].radius = radius[4];  // N3H1
				} else if(PDB[atomi].atomname[2] == 'E') {
					PDB[atomi].radius = radius[4];  // N3H1
				} else if((PDB[atomi].atomname[2] == 'D') && (!strncmp(PDB[atomi].res, "HIS", 3))) {
					PDB[atomi].radius = radius[4];  // N3H1
				} else if((PDB[atomi].atomname[2] == 'D') && (!strncmp(PDB[atomi].res, "ASN", 3))) {
					PDB[atomi].radius = radius[5];  // N3H2
				} else if((PDB[atomi].atomname[2] == 'E') && (!strncmp(PDB[atomi].res, "GLN", 3))) {
					PDB[atomi].radius = radius[5];  // N3H2
				} else if((PDB[atomi].atomname[2] == 'E') && (strncmp(PDB[atomi].res, "GLN", 3))) {
					PDB[atomi].radius = radius[4];  // N3H1
				} else if(PDB[atomi].atomname[2] == 'H') {
					PDB[atomi].radius = radius[5];  // N3H2
				} else if(PDB[atomi].atomname[2] == 'Z') {
					PDB[atomi].radius = radius[6];  // N4
				} else {
					PDB[atomi].radius = radius[5];  // N3H2,  N default
				}
			} else if(PDB[atomi].atomname[1] == 'C') {
				// CARBON
				if(PDB[atomi].atomname[2] == ' ') {
					PDB[atomi].radius = radius[0];  // C3H0 backbone carbon
				} else if((!strncmp(PDB[atomi].res, "ASP", 3)) && (PDB[atomi].atomname[2] == 'G')) {
					PDB[atomi].radius = radius[0];  // C3H0
				} else if((!strncmp(PDB[atomi].res, "GLU", 3)) && (PDB[atomi].atomname[2] == 'D')) {
					PDB[atomi].radius = radius[0];  // C3H0
				} else if((!strncmp(PDB[atomi].res, "ASN", 3)) && (PDB[atomi].atomname[2] == 'G')) {
					PDB[atomi].radius = radius[0];  // C3H0
				} else if((!strncmp(PDB[atomi].res, "GLN", 3)) && (PDB[atomi].atomname[2] == 'D')) {
					PDB[atomi].radius = radius[0];  // C3H0
				} else if((!strncmp(PDB[atomi].res, "ARG", 3)) && (PDB[atomi].atomname[2] == 'Z')) {
					PDB[atomi].radius = radius[0];  // C3H0
				} else if((!strncmp(PDB[atomi].res, "PHE", 3)) || (!strncmp(PDB[atomi].res, "HIS", 3))
			  || (!strncmp(PDB[atomi].res, "TYR", 3)) || (!strncmp(PDB[atomi].res, "TRP", 3))) {
					if((PDB[atomi].atomname[2] == 'A') || (PDB[atomi].atomname[2] == 'B')) {
						PDB[atomi].radius = radius[2];  // C4
					} else {
						PDB[atomi].radius = radius[1];  // C3H1
					}
				} else { // carbon is C4, aliphatic
					PDB[atomi].radius = radius[2];    // aliphatic carbon
				}
			} else {
				// default radius
				PDB[atomi].radius = radius[10];     // default for unknown atom;
			}
		}
	}
	return;
}

/******************************
 * subroutine index_protein
 * created 16/02/2001
 ******************************

// assigns all protein atoms to boxes within a cubic grid
// returns cube dimesions as number of boxes per side
// assigns values to global array box[]

 *****************************/

int index_protein(int PDBtot)
{
  int   xi;              // coordinate counter 0,1,2
  int   pdbi;            // atom counter 0 to (PDBtot-1)
  float cellsize = CELLSIZE;  // box dimensions
  float maxwidth;        // maximum protein dimension
  int   dim;             // dimension of DxDxD cube.
  int   boxi;            // box index
  int   startind;        // start point of box in PDBlist


  // ------ find global minimum and maximum -------
  for(xi=0; xi<3; ++xi) {
    globalmin[xi] = PDB[0].coor[xi];
    globalmax[xi] = PDB[0].coor[xi];
  }
  for (pdbi=0; pdbi<PDBtot; ++pdbi) {
    for (xi=0; xi<3; ++xi) {
      if (PDB[pdbi].coor[xi] < globalmin[xi]) {
	globalmin[xi] = PDB[pdbi].coor[xi];
      }
      if (PDB[pdbi].coor[xi] > globalmax[xi]) {
	globalmax[xi] = PDB[pdbi].coor[xi];
      }
    }
  }

  // ------ get largest dimension of protein -------
  maxwidth = 0.0;
  for(xi=0; xi<3; ++xi) {
    if((globalmax[xi]-globalmin[xi]) > maxwidth) {
      maxwidth = (globalmax[xi]-globalmin[xi]);
    }
  }
  // expand dimensions of cube by one (for round-off)
  dim = (int)(maxwidth/cellsize) + 1;

  // allocate memory for atomindex and PDBlist arrays
  box = calloc(dim*dim*dim, sizeof(struct atomindex));
  PDBlist = malloc(PDBtot*sizeof(int));
  if((!box) || (!PDBlist)) { 
    printf("memory allocation error for either 'box' or 'PDBlist'\n");
    exit(1);
  }

  // count entries per box, assign box number to atom
  for (pdbi=0; pdbi<PDBtot; ++pdbi) {
    boxi = (int)((PDB[pdbi].coor[0]-globalmin[0])/cellsize)*dim*dim
      + (int)((PDB[pdbi].coor[1]-globalmin[1])/cellsize)*dim
      + (int)((PDB[pdbi].coor[2]-globalmin[2])/cellsize);
    ++box[boxi].nument;
    PDB[pdbi].boxnum = boxi;
  }

  // assign start pointers for boxes in PDBlist
  startind = 0;
  for (boxi=0; boxi<dim*dim*dim; ++boxi) {
    box[boxi].first = startind;
    startind += box[boxi].nument;
  }

  // clear array (needed for recounting index)
  for(boxi=0; boxi<dim*dim*dim; ++boxi) {
    box[boxi].nument = 0;
  }

  // fill PDBlist index
  for (pdbi=0; pdbi<PDBtot; ++pdbi) {
    boxi = PDB[pdbi].boxnum;
    PDBlist[box[boxi].first +box[boxi].nument] = pdbi;
    ++box[boxi].nument;
  }
  return(dim);
}


/*********************************
 * subroutine get_contlist4
 *********************************/

// uses box index of protein atoms to find contacts in range of a given atom.
// requires global variable 'box[]'.
// checks previous atoms, keeps only those with non-zero contact area.

int get_contlist4(struct atom PDB[], int atomzero, struct contactlist contlist[], 
		  int PDBtot, float rado, int dim) 
{
  double sqrdist;             // distance squared between two points
  double neardist;            // max distance for contact between atom spheres
  int NC;                     // number of contacts
  int bai;                    // box atom counter
  int boxi;                   // current box number
  int atomj;                  // index number of atom in PDB
  int boxzero;                // box atomzero is in
  int i;                      // dummy counter for surrounding boxes
  long int currindex;

  NC = 0;
  // mark previously contacted atoms
  currindex = ca_index[atomzero];
  while(currindex != -1) {
    PDB[ca_rec[currindex].atom].done = 'C'; // makes contact
    currindex = ca_rec[currindex].prev;
  }

  // get pdb atom contacts from current and adjacent boxes
  boxzero = PDB[atomzero].boxnum;
  for(i=0; i<27; ++i) {
    // get up to 27 boxes surrounding current box
    boxi = boxzero +dim*dim*((i/9)-1) +dim*(((i/3)%3)-1) +(i%3)-1;
    if((boxi < 0) || (boxi >= dim*dim*dim)) continue;  // don't do boxes outside of range
    bai = 0;
    while(bai<box[boxi].nument) {
      atomj = PDBlist[box[boxi].first+bai]; 
      // check previous contacts
      if(PDB[atomj].done == 'Y') {
	++bai;
	continue;
      }
      sqrdist = (PDB[atomzero].coor[0]-PDB[atomj].coor[0])*(PDB[atomzero].coor[0]-PDB[atomj].coor[0]) 
	+ (PDB[atomzero].coor[1]-PDB[atomj].coor[1])*(PDB[atomzero].coor[1]-PDB[atomj].coor[1])
	+ (PDB[atomzero].coor[2]-PDB[atomj].coor[2])*(PDB[atomzero].coor[2]-PDB[atomj].coor[2]);
      neardist =  rado + PDB[atomj].radius +Rw;
      if((sqrdist < neardist*neardist) && (sqrdist != 0.0)) {
	// add atoms to list
	contlist[NC].index = atomj;
	contlist[NC].dist = sqrt(sqrdist);
	++NC;
      }
      ++bai;
    }
  }

  // reset atoms to 'done'
  currindex = ca_index[atomzero];
  while(currindex != -1) {
    PDB[ca_rec[currindex].atom].done = 'Y'; // re-mark as done
    currindex = ca_rec[currindex].prev;
  }
  return(NC);
}

/******************************
 * subroutine parse_commandline
 ******************************/

// parses the command line input for the program and
// sets options for command line switches.
// Sets global variables 'planedef' and 'nonbonded'.

void parse_commandline(int argc, char *argv[], char FILENAME[])
{
	int argi;  // argument counter
	showbonded = 'N';  // default
	planedef = 'R';    // default
	normalize = 'N';   // default
	argi = 1;

	if (argc < 2) {
		printf("\n==== Program  Vcontacts, version 1.2 ====\n");
		printf("This program calculates contact areas between atoms within a given structure file in pdb format.\n");
		printf("Atom-solvent as well as atom-atom contacts are calculated using a constrained Voronoi procedure. \n");
		printf("The contact dividing plane uses the radical plane of Gellatly and Finney as the default. Radii\n");
		printf("are read from the file 'radii.dat'.\n");
		printf("Usage: Vcontacts <pdb file>\n");
		printf("Optional switches:\n");
		printf("   -norm            normalizes contacts to a percent of non-bonded accessible surface\n");
		printf("   -planedef R      uses radical plane of Gellatly and Finney (DEFAULT)\n");
		printf("   -planedef X      uses extended radical plane of McConkey et al (matches solvent accessible surfaces)\n");
		printf("   -planedef B      uses the bisecting plane of the original Voronoi procedure\n");
		printf("\nPlease cite the following article as a reference:\n");
		printf("Quantification of protein surfaces, volumes and atom-atom contacts using a constrained Voronoi procedure\n");
		printf("B.J.McConkey, V. Sobolev, and M. Edelman (2002), Bioinformatics, 18:1365-1373.\n\n");
		exit(0);
	}

	while(argi < argc) {
		if(!strncmp(argv[argi], "-all", 4)) {
			showbonded = 'Y';
		} else if(!strncmp(argv[argi], "-norm", 5)) {
			normalize = 'Y';
		} else if(!strncmp(argv[argi], "-planedef",9)) {
			if((argv[argi][9] == 'X')||(argv[argi][9] == 'R')||(argv[argi][9] == 'B')) {
				planedef = argv[argi][9];
			} else if((argi+1) < argc) {
				if((argv[argi+1][0] == 'X')||(argv[argi+1][0] == 'R')||(argv[argi+1][0] == 'B')) {
					planedef = argv[argi+1][0];
					++argi;
				} else {
					printf("error: need to specify X, R, or B after -planedef switch\n");
					exit(1);
				}
			}	else {
				printf("error: need to specify X, R, or B after -planedef switch\n");
				exit(1);
			}
		} else {
			if(argv[argi][0] != '-') {
				// not a specified switch, so it's the filename
				strcpy(FILENAME, argv[argi]);
			}
		}
		++argi;
	}
	return;
}
