#ifndef DOMAIN_H
#define DOMAIN_H

/* Standard library includes - everything needed */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdbool.h>

/* ============================================================================
 * Geometry types (from geom.h)
 * ============================================================================ */
#ifndef PI
#define PI 3.14159265358979323846
#define twoPI (PI + PI)
#endif

typedef struct { float x, y, z; } Vec;
typedef struct { Vec A, B, C; } Mat;

/* Geometry function prototypes */
void  vinit(Vec*);
void  vcopy(Vec, Vec*);
void  vnorm(Vec*);
void  vave(Vec, Vec, Vec*);
void  vsum(Vec, Vec*);
void  vadd(Vec, Vec, Vec*);
void  vsub(Vec, Vec, Vec*);
void  vmul(Vec*, float);
void  vdiv(Vec*, float);
float vdif(Vec, Vec);
float vsqr(Vec);
float vmod(Vec);
float vdot(Vec, Vec);
void  vprod(Vec, Vec, Vec*);
void  VtoM(Vec, Vec, Vec, Mat*);
void  MmulV(Mat*, Vec, Vec*);
void  VmulM(Mat*, Vec, Vec*);

/* ============================================================================
 * PDB types (from pdbprot.h)
 * ============================================================================ */
typedef char Str4_[5];

typedef struct {
    int Atno;
    Str4_ Id;
    char Alt;
    char Aa;
    int Resno;
    char Rid;
    float X, Y, Z;
    float Occu, Bfact;
} Atom_;

typedef struct {
    int Don, Acc;
} Hbond_;

typedef enum { HELIX, SHEET, TURN } Sectype_;

typedef struct {
    Sectype_ Sectype;
    int No;
    Str4_ Id;
    int Beg, End;
    char Chid;
    char Begaa, Endaa;
    char Begrid, Endrid;
    int Type;
    int Strandno;
    Str4_ Thisat, Otherat;
    char Thisaa, Otheraa;
    int This, Other;
    char Thisrid, Otherid, Otherchid;
} Secstr_;

typedef struct {
    int No;
    int Pos1, Pos2;
    char Ch1, Ch2, Rid1, Rid2;
} Ssbond_;

typedef struct {
    Atom_ *Atoms;
    int Atomno;
    Secstr_ *Secs;
    int Secsno;
    Hbond_ *Hbonds;
    int Hbno;
    Ssbond_ *Ssbs;
    int Ssbno;
    int Aano;
    char Chid;
    char Type;
    char *Seq;
} Chain_;

typedef struct {
    char Header[41];
    char Date[10];
    Str4_ Pdbcode;
    char *Compound;
    char *Source;
    char Expdta[61];
    float Resol;
    Chain_ *Chains;
    int Chainno;
    Hbond_ *Hbonds;
    int Hbno;
    Ssbond_ *Ssbs;
    int Ssbno;
} Pdbentry_;

/* PDB function prototypes */
Pdbentry_* get_pdb(const char *Pdbfn, int Ca, int Strict);
void free_pdb(Pdbentry_ *Entry);

/* ============================================================================
 * Domain analysis configuration and types
 * ============================================================================ */
#define NCYCLES 200
#define NALLOC 2000
#define NACID 30
#define MINDOM 40
#define MAXDOM 400
#define DEFAULT_SPREAD 15
#define DEFAULT_NRUNS 1
#define IT 5

/* Structure definitions */
typedef struct {
    Vec v;
    float d;
    Vec cos;
} Tri;

typedef struct {
    int a, b;
    float c;
} Pairs;

typedef struct {
    char *res;
    float *acc;
    int *dom;
    int *rid;
    Vec *ca, *cb;
    int len;
} Seq;

typedef struct {
    int spread;
    int nruns;
    int subdom;
    bool sheet;
} DomainConfig;

/* Function prototypes */
bool domain_analysis_init(DomainConfig *config, int argc, char *argv[]);
void domain_analysis_cleanup(void);

int protin(Pdbentry_ *prot, Seq *seq, int id, Tri ***m, float z, int flip);
void add_cb(Seq *seq);
int copyca(Chain_ *pdb, Seq *s, int flip, float z);
void extend(Vec *res, int i, int j, int k, int new_idx);

void domain(float **mat, float **net, Vec *ca, float *ave, int n, 
            const DomainConfig *config);
void beta(float **net, float **mat, Vec *ca, float *ave, int n);
void ising(float **mat, float *ave, int n, int cycles, int limit, int smooth);
void smooth(float *dat, int n, int win);
void flatten(Vec *a, int n, int times);

int define(int is, int it, Seq *seq, Vec *cas, float *dom, 
           float **net, float **mat, int len, int *doms,
           const DomainConfig *config);

float compare(float **net, int *doms1, int *doms2, 
              int domn1, int domn2, int len, bool sheet);
float betacut(float **net, int *doms, int n);
int getbest(int **dd, int n);

void putpdb(Vec *atom, Seq *seq, FILE *out, int id, int subdom);
void loop(FILE *out, int id, int *n, Vec oldN, Vec oldC, 
          Vec cent, int offset);

void set_vect(Vec *a, Vec *b, Tri **m, int l);
void set_cbcb(Vec *a, Vec *b, Tri **m, int l);
void setframe(Vec a, Vec b, Vec c, Mat *frame);
void flipseq(Vec *ca, char *seq, float *acc, int n);

/* Memory management helpers */
float **alloc_matrix(int n);
void free_matrix(float **mat, int n);
int **alloc_int_matrix(int n);
void free_int_matrix(int **mat, int n);

#endif /* DOMAIN_H */
