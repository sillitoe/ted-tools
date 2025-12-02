/*
 * domain.c - Complete standalone protein domain identification tool
 * 
 * Self-contained implementation with all dependencies embedded.
 * Requires only standard C library and math library.
 *
 * Compilation:
 * cc -O2 -std=c99 -D_XOPEN_SOURCE=600 domain.c -o domain -lm
 *
 * Usage:
 * ./domain <pdb_file> [spread] [nruns] [subdom]
 */

#include "domain.h"

/* ============================================================================
 * GEOMETRY FUNCTIONS (from geom.c)
 * ============================================================================ */

void vinit(Vec *c) {
    c->x = 0.0f;
    c->y = 0.0f;
    c->z = 0.0f;
}

void vcopy(Vec b, Vec *c) {
    c->x = b.x;
    c->y = b.y;
    c->z = b.z;
}

void vnorm(Vec *c) {
    float d = sqrtf(c->x * c->x + c->y * c->y + c->z * c->z);
    if (d > 0.0f) {
        c->x /= d;
        c->y /= d;
        c->z /= d;
    }
}

void vave(Vec a, Vec b, Vec *c) {
    c->x = 0.5f * (a.x + b.x);
    c->y = 0.5f * (a.y + b.y);
    c->z = 0.5f * (a.z + b.z);
}

void vsum(Vec a, Vec *c) {
    c->x += a.x;
    c->y += a.y;
    c->z += a.z;
}

void vadd(Vec a, Vec b, Vec *c) {
    c->x = a.x + b.x;
    c->y = a.y + b.y;
    c->z = a.z + b.z;
}

void vsub(Vec a, Vec b, Vec *c) {
    c->x = a.x - b.x;
    c->y = a.y - b.y;
    c->z = a.z - b.z;
}

void vmul(Vec *c, float s) {
    c->x *= s;
    c->y *= s;
    c->z *= s;
}

void vdiv(Vec *c, float s) {
    c->x /= s;
    c->y /= s;
    c->z /= s;
}

float vdif(Vec a, Vec b) {
    Vec c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return sqrtf(c.x * c.x + c.y * c.y + c.z * c.z);
}

float vsqr(Vec a) {
    return (a.x * a.x + a.y * a.y + a.z * a.z);
}

float vmod(Vec a) {
    return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

float vdot(Vec a, Vec b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

void vprod(Vec a, Vec b, Vec *c) {
    c->x = a.y * b.z - b.y * a.z;
    c->y = a.z * b.x - b.z * a.x;
    c->z = a.x * b.y - b.x * a.y;
}

void VtoM(Vec a, Vec b, Vec c, Mat *M) {
    M->A.x = a.x; M->A.y = a.y; M->A.z = a.z;
    M->B.x = b.x; M->B.y = b.y; M->B.z = b.z;
    M->C.x = c.x; M->C.y = c.y; M->C.z = c.z;
}

void MmulV(Mat *m, Vec d, Vec *e) {
    Mat M = *m;
    Vec a = M.A, b = M.B, c = M.C;
    e->x = d.x * a.x + d.y * b.x + d.z * c.x;
    e->y = d.x * a.y + d.y * b.y + d.z * c.y;
    e->z = d.x * a.z + d.y * b.z + d.z * c.z;
}

void VmulM(Mat *m, Vec d, Vec *e) {
    Mat M = *m;
    Vec a = M.A, b = M.B, c = M.C;
    e->x = d.x * a.x + d.y * a.y + d.z * a.z;
    e->y = d.x * b.x + d.y * b.y + d.z * b.z;
    e->z = d.x * c.x + d.y * c.y + d.z * c.z;
}

/* ============================================================================
 * PDB I/O FUNCTIONS (from pdbprot.c)
 * ============================================================================ */

/* Convert 3-letter amino acid code to 1-letter */
static char aa_code31(const char *Aa3) {
    static char *Names3 = "ALAALBARGASNASPASXCYSCYHCSHCSSCYXGLNGLUGLXGLYHISILEILULEULYSMETPHEPROPR0PRZHYPSERTHRTRPTRYTYRVALUNKXXX";
    static char *Letters1 = "AARNDBCCCCCQEZGHIILKMFPPPPSTWWYVUX";
    char *Code;
    
    if (strlen(Aa3) < 3) return 'X';
    Code = strstr(Names3, Aa3);
    if (Code == NULL) return 'X';
    return Letters1[(Code - Names3) / 3];
}

/* Convert 1-letter code to 3-letter code */
char *aa_code13(char Aa1) {
    static char *Aacode[26] = {"ALA", "ASX", "CYS", "ASP", "GLU", "PHE",
        "GLY", "HIS", "ILE", "JJJ", "LYS", "LEU",
        "MET", "ASN", "OOO", "PRO", "GLN", "ARG",
        "SER", "THR", "UNK", "VAL", "TRP", "XXX",
        "TYR", "GLX"};
    
    if (Aa1 < 'A' || Aa1 > 'Z') return Aacode[23];
    return Aacode[Aa1 - 'A'];
}

/* Extract sequence from atom array */
static char *atom_seq(const Atom_ Atoms[], int Atomno, int Aano) {
    char *Seq = calloc(Aano + 1, sizeof(char));
    if (!Seq) return NULL;
    
    int j = 0;
    for (int i = 0; i < Atomno && j < Aano; j++) {
        Seq[j] = Atoms[i].Aa;
        char Rid = Atoms[i].Rid;
        int Resno = Atoms[i++].Resno;
        while (i < Atomno && Resno == Atoms[i].Resno && Rid == Atoms[i].Rid) i++;
    }
    return Seq;
}

/* Read PDB file - simplified version for C-alpha only */
Pdbentry_ *get_pdb(const char *Pdbfn, int Ca, int Strict) {
    FILE *Pdb = fopen(Pdbfn, "r");
    if (!Pdb) {
        fprintf(stderr, "Error: Cannot open %s\n", Pdbfn);
        return NULL;
    }
    
    Pdbentry_ *Entry = calloc(1, sizeof(Pdbentry_));
    if (!Entry) {
        fclose(Pdb);
        return NULL;
    }
    
    /* Initialize with safe defaults */
    Entry->Chains = NULL;
    Entry->Chainno = 0;
    Entry->Compound = calloc(61, 1);  /* Allocate proper space for compound name */
    Entry->Source = calloc(61, 1);
    strcpy(Entry->Compound, "Unknown");  /* Default compound name */
    strcpy(Entry->Expdta, "X-RAY DIFFRACTION");
    Entry->Resol = -1.0f;
    
    char Line[2000];
    Chain_ *Chains = NULL;
    Atom_ *Atoms = NULL;
    int Chainno = 0, Atomcnt = 0, Aano = 0;
    char Oldchain = '\0', Oldaa = '\0';
    int Oldresno = -9999;
    char Oldrid = ' ';
    int total_atoms_read = 0;  /* Track total atoms read for debugging */
    
    /* Read atom records */
    while (fgets(Line, sizeof(Line), Pdb)) {
        if (strncmp(Line, "HEADER", 6) == 0) {
            strncpy(Entry->Header, Line + 10, 40);
            Entry->Header[40] = '\0';
            if (strlen(Line) > 62) {
                strncpy(Entry->Pdbcode, Line + 62, 4);
                Entry->Pdbcode[4] = '\0';
            }
        }
        
        if (strncmp(Line, "COMPND", 6) == 0 && strlen(Line) > 10) {
            strncpy(Entry->Compound, Line + 10, 60);
            Entry->Compound[60] = '\0';
            /* Remove trailing newline if present */
            char *newline = strchr(Entry->Compound, '\n');
            if (newline) *newline = '\0';
        }
        
        if (strncmp(Line, "ATOM", 4) != 0) continue;
        
        total_atoms_read++;
        
        /* Check for chain change */
        char Chid = Line[21];
        if (Oldchain != '\0' && Oldchain != Chid) {
            /* Save previous chain */
            if (Atomcnt > 0) {  /* Only save if we have atoms */
                Chains = realloc(Chains, sizeof(Chain_) * (Chainno + 1));
                Chains[Chainno].Atoms = realloc(Atoms, Atomcnt * sizeof(Atom_));
                Chains[Chainno].Atomno = Atomcnt;
                Chains[Chainno].Aano = Aano;
                Chains[Chainno].Chid = Oldchain;
                Chains[Chainno].Type = 'A';
                Chains[Chainno].Seq = atom_seq(Chains[Chainno].Atoms, Atomcnt, Aano);
                Chains[Chainno].Secs = NULL;
                Chains[Chainno].Secsno = 0;
                Chains[Chainno].Hbonds = NULL;
                Chains[Chainno].Hbno = 0;
                Chains[Chainno].Ssbs = NULL;
                Chains[Chainno].Ssbno = 0;
                Chainno++;
            }
            
            /* Start new chain */
            Atoms = NULL;
            Atomcnt = 0;
            Aano = 0;
            Oldresno = -9999;
            Oldrid = ' ';
        }
        Oldchain = Chid;
        
        /* Parse atom - skip if not CA and Ca flag is set */
        char Atname[5];
        if (Line[12] == ' ') {
            strncpy(Atname, Line + 13, 3);
            Atname[3] = '\0';
        } else {
            strncpy(Atname, Line + 12, 4);
            Atname[4] = '\0';
        }
        
        /* Remove trailing spaces */
        for (int i = strlen(Atname) - 1; i >= 0 && Atname[i] == ' '; i--) {
            Atname[i] = '\0';
        }
        
        if (Ca && strcmp(Atname, "CA") != 0) continue;
        
        /* Allocate space */
        if (Atomcnt % 256 == 0) {
            Atoms = realloc(Atoms, (Atomcnt + 256) * sizeof(Atom_));
        }
        
        /* Parse atom data */
        sscanf(Line, "ATOM%d", &Atoms[Atomcnt].Atno);
        strcpy(Atoms[Atomcnt].Id, Atname);
        Atoms[Atomcnt].Alt = Line[16];
        
        char Aa3[4];
        strncpy(Aa3, Line + 17, 3);
        Aa3[3] = '\0';
        Atoms[Atomcnt].Aa = aa_code31(Aa3);
        
        if (Strict && Atoms[Atomcnt].Aa == 'X') continue;
        
        sscanf(Line + 22, "%4d", &Atoms[Atomcnt].Resno);
        Atoms[Atomcnt].Rid = Line[26];
        
        sscanf(Line + 30, "%8f%8f%8f%6f%6f",
               &Atoms[Atomcnt].X, &Atoms[Atomcnt].Y, &Atoms[Atomcnt].Z,
               &Atoms[Atomcnt].Occu, &Atoms[Atomcnt].Bfact);
        
        /* Track residue count */
        if (Atoms[Atomcnt].Resno != Oldresno || 
            Atoms[Atomcnt].Rid != Oldrid ||
            Atoms[Atomcnt].Aa != Oldaa) {
            Oldresno = Atoms[Atomcnt].Resno;
            Oldrid = Atoms[Atomcnt].Rid;
            Oldaa = Atoms[Atomcnt].Aa;
            Aano++;
        }
        
        Atomcnt++;
    }
    
    /* Save last chain */
    if (Atomcnt > 0) {
        Chains = realloc(Chains, sizeof(Chain_) * (Chainno + 1));
        Chains[Chainno].Atoms = realloc(Atoms, Atomcnt * sizeof(Atom_));
        Chains[Chainno].Atomno = Atomcnt;
        Chains[Chainno].Aano = Aano;
        Chains[Chainno].Chid = Oldchain;
        Chains[Chainno].Type = 'A';
        Chains[Chainno].Seq = atom_seq(Chains[Chainno].Atoms, Atomcnt, Aano);
        Chains[Chainno].Secs = NULL;
        Chains[Chainno].Secsno = 0;
        Chains[Chainno].Hbonds = NULL;
        Chains[Chainno].Hbno = 0;
        Chains[Chainno].Ssbs = NULL;
        Chains[Chainno].Ssbno = 0;
        Chainno++;
    }
    
    Entry->Chains = Chains;
    Entry->Chainno = Chainno;
    Entry->Hbonds = NULL;
    Entry->Hbno = 0;
    Entry->Ssbs = NULL;
    Entry->Ssbno = 0;
    
    fclose(Pdb);
    
    /* Debug output */
    fprintf(stderr, "Debug: Read %d total ATOM lines\n", total_atoms_read);
    fprintf(stderr, "Debug: Found %d chains\n", Entry->Chainno);
    if (Entry->Chainno > 0) {
        fprintf(stderr, "Debug: First chain has %d CA atoms\n", Entry->Chains[0].Atomno);
    }
    
    return Entry;
}

/* Free PDB entry */
void free_pdb(Pdbentry_ *Entry) {
    if (!Entry) return;
    
    for (int i = 0; i < Entry->Chainno; i++) {
        if (Entry->Chains[i].Seq) free(Entry->Chains[i].Seq);
        if (Entry->Chains[i].Atoms) free(Entry->Chains[i].Atoms);
        if (Entry->Chains[i].Secs) free(Entry->Chains[i].Secs);
        if (Entry->Chains[i].Hbonds) free(Entry->Chains[i].Hbonds);
        if (Entry->Chains[i].Ssbs) free(Entry->Chains[i].Ssbs);
    }
    if (Entry->Chains) free(Entry->Chains);
    if (Entry->Ssbs) free(Entry->Ssbs);
    if (Entry->Hbonds) free(Entry->Hbonds);
    free(Entry->Compound);
    free(Entry->Source);
    free(Entry);
}

/* ============================================================================
 * MEMORY MANAGEMENT HELPERS
 * ============================================================================ */

float **alloc_matrix(int n) {
    float **mat = calloc(n, sizeof(float*));
    if (!mat) return NULL;
    
    for (int i = 0; i < n; i++) {
        mat[i] = calloc(n, sizeof(float));
        if (!mat[i]) {
            for (int j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void free_matrix(float **mat, int n) {
    if (!mat) return;
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

int **alloc_int_matrix(int n) {
    int **mat = calloc(n, sizeof(int*));
    if (!mat) return NULL;
    
    for (int i = 0; i < n; i++) {
        mat[i] = calloc(n, sizeof(int));
        if (!mat[i]) {
            for (int j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void free_int_matrix(int **mat, int n) {
    if (!mat) return;
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

/* ============================================================================
 * DOMAIN ANALYSIS ALGORITHM
 * ============================================================================ */

/* Configuration initialization */
bool domain_analysis_init(DomainConfig *config, int argc, char *argv[]) {
    config->spread = DEFAULT_SPREAD;
    config->nruns = DEFAULT_NRUNS;
    config->subdom = 0;
    config->sheet = false;
    
    if (argc > 2) {
        if (sscanf(argv[2], "%d", &config->spread) != 1) {
            fprintf(stderr, "Error: Invalid spread parameter\n");
            return false;
        }
    }
    if (argc > 3) {
        if (sscanf(argv[3], "%d", &config->nruns) != 1) {
            fprintf(stderr, "Error: Invalid nruns parameter\n");
            return false;
        }
    }
    if (argc > 4) {
        if (sscanf(argv[4], "%d", &config->subdom) != 1) {
            fprintf(stderr, "Error: Invalid subdom parameter\n");
            return false;
        }
    }
    
    printf("Configuration:\n");
    printf("  spread = %d\n  nruns = %d\n  subdom = %d\n\n", 
           config->spread, config->nruns, config->subdom);
    
    return true;
}

/* Check for chain breaks */
static void check_chain_breaks(const Seq *seq) {
    for (int i = 1; i < seq->len; i++) {
        float d = vdif(seq->ca[i], seq->ca[i+1]);
        if (d > 8.0f) {
            printf("*NB* chain BREAK at %d/%d (distance: %.2f)\n", i, i+1, d);
        }
    }
}

/* Find best matching pair in overlap matrix */
int getbest(int **dd, int n) {
    int max = 0;
    int row = 0, col = 0;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dd[i][j] > max) {
                max = dd[i][j];
                row = j;
                col = i;
            }
        }
    }
    
    if (max == 0) return 0;
    
    for (int i = 0; i < n; i++) {
        dd[col][i] = -1;
        dd[i][row] = -1;
    }
    
    return max;
}

/* Calculate beta sheet splitting error */
float betacut(float **net, int *doms, int n) {
    int in[200] = {0};
    float cut = 0.0f;
    
    /* Check for domain indices that exceed array bounds */
    for (int i = 1; i <= n; i++) {
        if (doms[i] >= 0 && doms[i] < 200) {
            in[doms[i]]++;
        }
    }
    
    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j <= n; j++) {
            if (net[i][j] < 0.0f) continue;
            if (doms[i] == doms[j]) continue;
            if (doms[i] >= 200 || doms[j] >= 200) continue;  /* Bounds check */
            if (doms[i] < 0 || doms[j] < 0) continue;        /* Negative check */
            if (in[doms[j]] <= MINDOM || in[doms[i]] <= MINDOM) continue;
            cut += net[i][j];
        }
    }
    
    return cut;
}

/* Compare two domain assignments */
float compare(float **net, int *doms1, int *doms2, 
              int domn1, int domn2, int len, bool sheet) {
    int maxdd = 0;
    for (int i = 1; i <= len; i++) {
        if (doms1[i] > maxdd) maxdd = doms1[i];
        if (doms2[i] > maxdd) maxdd = doms2[i];
    }
    maxdd++;  /* Need space for indices 0 to maxdd */
    
    int **dd = alloc_int_matrix(maxdd);
    if (!dd) return -1.0f;
    
    /* Initialize matrix to zero */
    for (int i = 0; i < maxdd; i++) {
        for (int j = 0; j < maxdd; j++) {
            dd[i][j] = 0;
        }
    }
    
    /* Build overlap matrix */
    for (int i = 1; i <= len; i++) {
        if (doms1[i] >= 0 && doms1[i] < maxdd && 
            doms2[i] >= 0 && doms2[i] < maxdd) {
            dd[doms1[i]][doms2[i]]++;
        }
    }
    
    int match = 0;
    int best;
    while ((best = getbest(dd, maxdd)) > 0) {
        match += best;
    }
    
    float pct = 100.0f * (float)match / (float)len;
    
    if (sheet) {
        float bad = sqrtf((float)len);
        float cuts1 = betacut(net, doms1, len);
        printf("Error in beta split 1 = %7.2f\n", cuts1);
        float cuts2 = betacut(net, doms2, len);
        printf("Error in beta split 2 = %7.2f\n", cuts2);
        printf("(error limit = %7.3f)\n", bad);
        if (cuts1 > bad || cuts2 > bad) {
            pct = -pct;
        }
    }
    
    printf("%d on %d domains, %3d pct. agree\n", domn1, domn2, (int)pct);
    
    free_int_matrix(dd, maxdd);
    return pct;
}

/* Ising model optimization */
void ising(float **mat, float *ave, int n, int cycles, int limit, int smooth) {
    float *old = calloc((n + 1), sizeof(float));
    float *new = calloc((n + 1), sizeof(float));
    int **list = calloc((n + 1), sizeof(int*));
    
    if (!old || !new || !list) {
        free(old);
        free(new);
        free(list);
        return;
    }
    
    for (int j = 1; j <= n; j++) {
        old[j] = new[j] = ave[j];
    }
    
    for (int i = 1; i <= n; i++) {
        int count = 0;
        for (int j = 1; j <= n; j++) {
            if (i != j && mat[i][j] >= 0.0f) count++;
        }
        
        list[i] = calloc((count + 1), sizeof(int));
        if (!list[i]) {
            for (int k = 1; k < i; k++) free(list[k]);
            free(old);
            free(new);
            free(list);
            return;
        }
        
        int k = 1;
        for (int j = 1; j <= n; j++) {
            if (i != j && mat[i][j] >= 0.0f) {
                list[i][k++] = j;
            }
        }
        list[i][0] = k;
    }
    
    float shift = 1.0f;
    float step = 2.0f / (float)n;
    
    int k;
    for (k = 1; k < cycles; k++) {
        for (int j = 1; j <= n; j++) {
            float sum = 0.0f;
            float mean = 0.0f;
            float wt = 0.0f;
            int in = 0;
            
            for (int ii = 1; ii < list[j][0]; ii++) {
                int i = list[j][ii];
                if (limit && abs(i - j) > limit) continue;
                
                if (old[i] > old[j]) sum += mat[i][j];
                if (old[i] < old[j]) sum -= mat[i][j];
                
                if (old[i] > 0.0f) {
                    mean += old[i] * mat[i][j];
                    wt += mat[i][j];
                    in++;
                }
            }
            
            if (sum > 0.0f) new[j] += shift;
            if (sum < 0.0f) new[j] -= shift;
            if (in && (old[j] < 0.0f || smooth)) {
                new[j] = mean / wt;
            }
        }
        
        float rms = 0.0f;
        
        for (int j = 1; j <= n; j++) {
            float a = ave[j];
            ave[j] = 0.5f * (new[j] + old[j]);
            a -= ave[j];
            rms += a * a;
            old[j] = new[j];
        }
        
        rms = sqrtf(rms / (float)n);
        
        if (rms < 0.000001f || k > n / 2) {
            shift -= step;
        }
        if (shift < 0.0f) break;
    }
    
    printf("%d iterations\n", k);
    
    free(old);
    free(new);
    for (int i = 1; i <= n; i++) {
        free(list[i]);
    }
    free(list);
}

/* Smooth array using median filter */
void smooth(float *dat, int n, int win) {
    int w = win / 2;
    /* compute exact buffer size: code writes up to index n + 2*win + 1 */
    int buf_size = n + win * 2 + 2; /* allows indices 0 .. n+2*win+1 */
    float *old = calloc(buf_size, sizeof(float));
    float *new = calloc(buf_size, sizeof(float));
    /* p needs to be at least win+1 because code accesses p[w] */
    int *p = calloc((win + 1), sizeof(int));
    
    if (!old || !new || !p) {
        free(old);
        free(new);
        free(p);
        return;
    }
    
    for (int i = 1; i <= n; i++) {
        old[i + win] = dat[i];
    }
    for (int i = 0; i <= win; i++) {
        old[i] = new[i] = dat[win - i + 1];
        old[i + n + win + 1] = new[i + n + win + 1] = dat[n - i];
    }
    
    for (int k = 1; k < win; k++) {
        float rms = 0.0f;
        
        for (int i = w; i <= n + win + w; i++) {
            for (int j = 0; j < win; j++) {
                p[j] = j;
            }
            for (int j = 0; j < win - 1; j++) {
                for (int m = j + 1; m < win; m++) {
                    if (old[i - w + p[j]] > old[i - w + p[m]]) {
                        int tmp = p[j];
                        p[j] = p[m];
                        p[m] = tmp;
                    }
                }
            }
            new[i] = old[i + p[w] - w];
        }
        
        dat[0] = new[win];
        for (int i = 1; i <= n; i++) {
            float a = old[i + win] - new[i + win];
            rms += a * a;
            dat[i] = new[i + win];
        }
        
        rms = sqrtf(rms / (float)n);
        if (rms < 0.000001f) break;
        
        /* copy full valid range (0 .. n+2*win+1) */
        for (int i = 0; i <= n + win * 2 + 1; i++) {
            old[i] = new[i];
        }
    }
    
    free(old);
    free(new);
    free(p);
}

/* Flatten/smooth coordinate trajectory */
void flatten(Vec *a, int n, int times) {
    for (int i = 0; i < times; i++) {
        Vec p1, p2, p3;
        vcopy(a[1], &p1);
        vcopy(a[2], &p2);
        vcopy(a[3], &p3);
        
        for (int j = 2; j < n; j++) {
            Vec q;
            vinit(&q);
            vsum(p1, &q);
            vsum(p2, &q);
            vsum(p3, &q);
            vdiv(&q, 3.0f);
            
            vcopy(a[j], &p1);
            vcopy(a[j + 1], &p2);
            vcopy(a[j + 2], &p3);
            vcopy(q, a + j);
        }
    }
}

/* Beta sheet detection */
void beta(float **net, float **mat, Vec *ca, float *ave, int n) {
    const float cut = 7.5f;
    
    printf("\n");
    
    for (int i = 1; i <= n; i++) {
        ave[i] = (float)i;
    }
    
    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j <= n; j++) {
            mat[i][j] = mat[j][i] = vdif(ca[i], ca[j]);
        }
    }
    
    for (int i = 0; i < n + 2; i++) {
        for (int j = 0; j < n + 2; j++) {
            net[i][j] = 0.0f;
        }
    }
    
    for (int i = 2; i < n; i++) {
        int p[3] = {0, 0, 0};
        int q[3] = {0, 0, 0};
        int s[3] = {i - 1, i, i + 1};
        float dmin;
        
        dmin = 999.9f;
        for (int j = 2; j < n; j++) {
            float a = mat[i][j];
            if (abs(i - j) < 4 || a > cut || a > dmin) continue;
            p[1] = j;
            dmin = a;
        }
        
        dmin = 999.9f;
        for (int j = 2; j < n; j++) {
            float a = mat[i][j];
            if (abs(i - j) < 4 || a > cut || a > dmin || abs(j - p[1]) < 6) continue;
            q[1] = j;
            dmin = a;
        }
        
        if (!p[1] || !q[1]) continue;
        
        if (fminf(mat[s[0]][p[1] - 1], mat[s[2]][p[1] + 1]) <
            fminf(mat[s[2]][p[1] - 1], mat[s[0]][p[1] + 1])) {
            p[0] = p[1] - 1;
            p[2] = p[1] + 1;
        } else {
            p[0] = p[1] + 1;
            p[2] = p[1] - 1;
        }
        
        if (fminf(mat[s[0]][q[1] - 1], mat[s[2]][q[1] + 1]) <
            fminf(mat[s[2]][q[1] - 1], mat[s[0]][q[1] + 1])) {
            q[0] = q[1] - 1;
            q[2] = q[1] + 1;
        } else {
            q[0] = q[1] + 1;
            q[2] = q[1] - 1;
        }
        
        if (mat[s[0]][p[0]] > cut || mat[s[2]][p[2]] > cut ||
            mat[s[0]][q[0]] > cut || mat[s[2]][q[2]] > cut) continue;
        
        bool valid = true;
        for (int k = 0; k < 3 && valid; k++) {
            for (int m = 0; m < 3 && valid; m++) {
                if (abs(s[k] - p[m]) < 4 || abs(s[k] - q[m]) < 4 ||
                    abs(q[k] - p[m]) < 4) {
                    valid = false;
                }
            }
        }
        if (!valid) continue;
        
        printf("sheet: %4d %4d %4d\n", p[1], s[1], q[1]);
        
        for (int k = 0; k < 3; k++) {
            net[s[k]][q[k]] += 1.0f;
            net[q[k]][s[k]] += 1.0f;
            net[s[k]][p[k]] += 1.0f;
            net[p[k]][s[k]] += 1.0f;
        }
        
        net[s[1]][s[0]] += 1.0f;
        net[s[0]][s[1]] += 1.0f;
        net[s[1]][s[2]] += 1.0f;
        net[s[2]][s[1]] += 1.0f;
        net[p[1]][p[0]] += 1.0f;
        net[p[0]][p[1]] += 1.0f;
        net[p[1]][p[2]] += 1.0f;
        net[p[2]][p[1]] += 1.0f;
        net[q[1]][q[0]] += 1.0f;
        net[q[0]][q[1]] += 1.0f;
        net[q[1]][q[2]] += 1.0f;
        net[q[2]][q[1]] += 1.0f;
    }
    
    printf("\n");
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (net[i][j] > 0.5f) {
                net[i][i] += 1.0f;
            } else if (i != j) {
                net[i][j] = -1.0f;
            }
        }
    }
}

/* Calculate domain boundaries */
void domain(float **mat, float **net, Vec *ca, float *ave, int n,
            const DomainConfig *config) {
    if (config->sheet) {
        printf("setting sheet\n");
        ising(net, ave, n, n, 0, 1);
    }
    
    printf("defining domains\n");
    
    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j <= n; j++) {
            mat[i][j] = mat[j][i] = vdif(ca[i], ca[j]);
        }
    }
    
    for (int i = 1; i <= n; i++) {
        mat[i][i] = 0.0f;
        for (int j = 1; j <= n; j++) {
            if (i == j) continue;
            
            float d = mat[i][j];
            mat[i][j] = (float)config->spread / d;
            
            if (d > (float)config->spread) {
                mat[i][j] = -1.0f;
            } else {
                mat[i][i] += 1.0f;
            }
            
            if (ave[i] < 0.0f && ave[j] < 0.0f) {
                mat[i][j] = -1.0f;
            }
        }
    }
    
    ising(mat, ave, n, n, 0, 0);
    smooth(ave, n, 21);
}

/* PDB output functions */
void putpdb(Vec *atom, Seq *seq, FILE *out, int id, int subdom);
void loop(FILE *out, int id, int *n, Vec oldN, Vec oldC, Vec cent, int offset);

/* Define domain boundaries */
int define(int is, int it, Seq *seq, Vec *cas, float *dom,
           float **net, float **mat, int len, int *doms,
           const DomainConfig *config) {
    printf("\n\n**** CYCLE %d ****\n\n", it);
    
    float *tmp = calloc((len + 2), sizeof(float));
    int *p = calloc((len + 1), sizeof(int));
    
    if (!tmp || !p) {
        free(tmp);
        free(p);
        return 0;
    }
    
    domain(mat, net, cas, dom, len, config);
    
    for (int i = 1; i <= len; i++) {
        dom[i] = -dom[i];
    }
    dom[0] = 999.9f;
    
    for (int i = 0; i <= len; i++) {
        p[i] = i;
    }
    for (int i = 1; i <= len; i++) {
        for (int j = i + 1; j <= len; j++) {
            if (dom[p[i]] < dom[p[j]]) {
                int temp = p[i];
                p[i] = p[j];
                p[j] = temp;
            }
        }
    }
    
    seq->dom[p[0]] = 0;
    int n = 0;
    int in[200] = {0};
    
    for (int i = 1; i <= len; i++) {
        float d1 = dom[p[i - 1]];
        float d2 = dom[p[i]];
        
        if (d1 - d2 > 1.5f) {
            n++;
            in[n] = 0;
        }
        seq->dom[p[i]] = n;
        in[n]++;
    }
    
    int ndom = n;
    float f = (float)(n - 1);
    for (int i = 1; i <= len; i++) {
        seq->acc[i] = (n > 1) ? (float)(seq->dom[i] - 1) / f : 0.0f;
    }
    
    printf("\n");
    
    int redo = 0;
    if (n == 1) {
        printf("assign   1 domain\n");
        *doms = 1;
        free(tmp);
        free(p);
        return 0;
    }
    
    printf("%d domains\n", ndom);
    *doms = 0;
    
    for (int i = 1; i <= n; i++) {
        printf("   %5d = %d res.\n", i, in[i]);
        if (in[i] < MINDOM) {
            redo = 1;
        } else {
            (*doms)++;
        }
    }
    
    if (*doms == 0) {
        printf("NO domains\n");
        free(tmp);
        free(p);
        return 0;
    }
    
    int last = (!redo || it == IT);
    int mark = (last && !is);
    
    if (mark && !config->nruns) {
        printf("assign   %d domains\n", *doms);
        
        char filename[255];
        if (config->subdom) {
            sprintf(filename, "dom%d-0.out", config->subdom);
        } else {
            sprintf(filename, "dom0.out");
        }
        
        FILE *out = fopen(filename, "w");
        if (out) {
            putpdb(seq->ca, seq, out, 0, config->subdom);
            fclose(out);
        }
        
        int j = 0;
        for (int i = 1; i <= n; i++) {
            if (in[i] < MINDOM) continue;
            
            printf("assign   %5d = %d res.\n", i, in[i]);
            j++;
            
            if (config->subdom) {
                sprintf(filename, "dom%d-%d.out", config->subdom, j);
            } else {
                sprintf(filename, "dom%d.out", j);
            }
            
            out = fopen(filename, "w");
            if (out) {
                putpdb(seq->ca, seq, out, i, config->subdom);
                fclose(out);
            }
        }
    }
    
    printf("\n");
    n = 1;
    if (mark && !config->nruns) printf("assign   ");
    printf("segment 1 in domain %d = 1 (%c%d)", 
           seq->dom[1], seq->res[1], seq->rid[1]);
    
    for (int i = 2; i <= len; i++) {
        if (seq->dom[i - 1] == seq->dom[i]) continue;
        
        n++;
        printf(" --- %d (%c%d) \n", i - 1, seq->res[i - 1], seq->rid[i - 1]);
        if (mark && !config->nruns) printf("assign   ");
        printf("segment %d in domain %d = %d (%c%d)", 
               n, seq->dom[i], i, seq->res[i], seq->rid[i]);
    }
    
    printf(" --- %d (%c%d)\n", len, seq->res[len], seq->rid[len]);
    if (mark && !config->nruns) printf("assign   ");
    printf("%d segments\n", n);
    
    for (int i = 1; i <= len; i++) {
        if (in[seq->dom[i]] < MINDOM) {
            dom[i] = -1.0f;
        } else {
            dom[i] = 100.0f * (float)(seq->dom[i]);
        }
    }
    
    free(tmp);
    free(p);
    return redo;
}

/* PDB output */
void putpdb(Vec *atom, Seq *seq, FILE *out, int id, int subdom) {
    int n = 0;
    int len = seq->len;
    int insert = seq->rid[len] + 100;
    int offset = subdom ? 100 + insert : insert;
    int Ndom = 0, Cdom = 0;
    int start = 0, end = len;
    
    char aa3[80];
    strcpy(aa3,
        "ALAASXCYSASPGLUPHEGLYHISILEACELYSLEUMETASNPCAPROGLNARGSERTHRUNKVALTRPXXXTYRGLX");
    
    for (int i = len; i > 0; i--) {
        if (id && id == seq->dom[i] && seq->rid[i] < insert) {
            end = i;
            break;
        }
    }
    
    for (int i = 1; i <= end; i++) {
        char aaa[4];
        strncpy(aaa, aa3 + 3 * (seq->res[i] - 'A'), 3);
        aaa[3] = '\0';
        
        if (id && id != seq->dom[i]) {
            Cdom = i;
            continue;
        }
        
        if (!start && seq->rid[i] >= insert) continue;
        start = 1;
        
        if (Ndom && Cdom) {
            Vec cent;
            vinit(&cent);
            for (int j = Ndom + 1; j < Cdom; j++) {
                vsum(atom[j], &cent);
            }
            vdiv(&cent, (float)(Cdom - Ndom - 1));
            loop(out, id, &n, atom[Ndom], atom[Cdom], cent, offset);
            Ndom = Cdom = 0;
        }
        
        n++;
        fprintf(out, "ATOM%7d  CA  %s %c%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n",
                n, aaa, 'A' + id, seq->rid[i], 
                atom[i].x, atom[i].y, atom[i].z, seq->acc[i]);
        
        if (id && seq->dom[i] != seq->dom[i + 1]) {
            Ndom = i + 1;
        }
    }
    
    fprintf(out, "TER\n");
}

void loop(FILE *out, int id, int *n, Vec oldN, Vec oldC, Vec cent, int offset) {
    *n = *n + 1;
    fprintf(out, "ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
            *n, 'A' + id, *n + offset, oldN.x, oldN.y, oldN.z);
    
    if (vdif(oldN, oldC) > 5.0f) {
        Vec newN, newC, ave;
        
        vsub(cent, oldN, &newN);
        vsub(cent, oldC, &newC);
        vnorm(&newN);
        vnorm(&newC);
        vmul(&newN, 3.8f);
        vmul(&newC, 3.8f);
        vadd(oldN, newN, &newN);
        vadd(oldC, newC, &newC);
        vave(newN, newC, &ave);
        
        if (vdif(newN, newC) < 3.0f) {
            *n = *n + 1;
            fprintf(out, "ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
                    *n, 'A' + id, *n + offset, ave.x, ave.y, ave.z);
        } else {
            vave(ave, cent, &cent);
            loop(out, id, n, newN, newC, cent, offset);
        }
    }
    
    *n = *n + 1;
    fprintf(out, "ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
            *n, 'A' + id, *n + offset, oldC.x, oldC.y, oldC.z);
}

/* Protein structure initialization */
void extend(Vec *res, int i, int j, int k, int new_idx);
void add_cb(Seq *seq);
int copyca(Chain_ *pdb, Seq *s, int flip, float z);
void set_vect(Vec *a, Vec *b, Tri **m, int l);
void set_cbcb(Vec *a, Vec *b, Tri **m, int l);
void setframe(Vec a, Vec b, Vec c, Mat *frame);
void flipseq(Vec *ca, char *seq, float *acc, int n);

int protin(Pdbentry_ *prot, Seq *seq, int id, Tri ***m, float z, int flip) {
    (void)id;  /* Unused parameter */
    
    /* Validate input */
    if (!prot || !seq) {
        fprintf(stderr, "Error: NULL protein or sequence structure\n");
        return 0;
    }
    
    if (prot->Chainno == 0 || !prot->Chains) {
        fprintf(stderr, "Error: No chains found in protein structure\n");
        return 0;
    }
    
    if (prot->Chains[0].Aano == 0 || !prot->Chains[0].Atoms) {
        fprintf(stderr, "Error: No atoms found in first chain\n");
        return 0;
    }
    
    int len = copyca(prot->Chains, seq, flip, z);
    
    if (len < 3) {
        fprintf(stderr, "Error: Protein too short (%d residues), minimum 3 required\n", len);
        return 0;
    }
    
    Vec cog;
    vinit(&cog);
    for (int i = 1; i <= len; i++) {
        vsum(seq->ca[i], &cog);
    }
    vdiv(&cog, (float)len);
    
    for (int i = 1; i <= len; i++) {
        Vec temp;
        vsub(seq->ca[i], cog, &temp);
        vcopy(temp, seq->ca + i);
    }
    
    Tri **mat = calloc((len + 2), sizeof(Tri*));
    if (!mat) return 0;

    for (int i = 0; i <= len + 1; i++) {
        mat[i] = calloc((len + 2), sizeof(Tri));
        if (!mat[i]) {
            for (int j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return 0;
        }
    }
    
    add_cb(seq);
    set_vect(seq->ca, seq->cb, mat, len);
    set_cbcb(seq->ca, seq->cb, mat, len);
    
    *m = mat;
    return len;
}

void add_cb(Seq *seq) {
    for (int i = 1; i <= seq->len; i++) {
        Vec n, c, b;
        float bond = 3.0f;
        
        vsub(seq->ca[i], seq->ca[i - 1], &n);
        vnorm(&n);
        vsub(seq->ca[i], seq->ca[i + 1], &c);
        vnorm(&c);
        vadd(n, c, &b);
        vnorm(&b);
        vmul(&b, bond);
        vadd(seq->ca[i], b, &seq->cb[i]);
    }
}

int copyca(Chain_ *pdb, Seq *s, int flip, float z) {
    int n = pdb->Aano;
    
    /* Validate minimum chain length */
    if (n < 3) {
        fprintf(stderr, "Error: Chain too short (%d residues)\n", n);
        return 0;
    }
    
    char *seq = calloc((n + 3), sizeof(char));
    float *acc = calloc((n + 3), sizeof(float));
    int *dom = calloc((n + 3), sizeof(int));
    int *rid = calloc((n + 3), sizeof(int));
    Vec *ca = calloc((n + 3), sizeof(Vec));
    Vec *cb = calloc((n + 3), sizeof(Vec));
    
    if (!seq || !acc || !dom || !rid || !ca || !cb) {
        free(seq);
        free(acc);
        free(dom);
        free(rid);
        free(ca);
        free(cb);
        return 0;
    }
    
    for (int i = 0; i < n; i++) {
        ca[i + 1].x = pdb->Atoms[i].X;
        ca[i + 1].y = pdb->Atoms[i].Y;
        ca[i + 1].z = pdb->Atoms[i].Z;
        acc[i + 1] = pdb->Atoms[i].Bfact;
        rid[i + 1] = pdb->Atoms[i].Resno;
        seq[i + 1] = pdb->Atoms[i].Aa;
        
        if (seq[i + 1] < 'A' || seq[i + 1] > 'Z') {
            printf("*NB* funny aa = %c\n", seq[i + 1]);
            seq[i + 1] = 'X';
        }
    }
    
    seq[0] = 'n';
    extend(ca, 3, 2, 1, 0);
    extend(ca, n - 2, n - 1, n, n + 1);
    seq[n + 1] = 'c';
    seq[n + 2] = '\0';
    
    for (int i = 0; i <= n + 1; i++) {
        ca[i].z *= z;
    }
    
    if (flip) {
        flipseq(ca, seq, acc, n);
    }
    
    s->res = seq;
    s->acc = acc;
    s->dom = dom;
    s->rid = rid;
    s->ca = ca;
    s->cb = cb;
    s->len = n;
    
    return n;
}

void extend(Vec *res, int i, int j, int k, int new_idx) {
    Vec m, v;
    vave(res[j], res[k], &m);
    vsub(m, res[i], &v);
    vadd(m, v, &res[new_idx]);
}

void flipseq(Vec *ca, char *seq, float *acc, int n) {
    for (int i = 0; i <= n / 2; i++) {
        int j = n + 1 - i;
        
        Vec r = ca[i];
        ca[i] = ca[j];
        ca[j] = r;
        
        char c = seq[i];
        seq[i] = seq[j];
        seq[j] = c;
        
        float a = acc[i];
        acc[i] = acc[j];
        acc[j] = a;
    }
}

void set_vect(Vec *a, Vec *b, Tri **m, int l) {
    (void)b;  /* Unused parameter */
    for (int i = 1; i <= l; i++) {
        Mat frame;
        setframe(a[i - 1], a[i], a[i + 1], &frame);
        
        for (int j = 1; j <= l; j++) {
            m[i][j].d = vdif(a[i], a[j]);
            vinit(&(m[i][j].v));
            
            if (i == j) continue;
            
            Vec s;
            vsub(a[j], a[i], &s);
            VmulM(&frame, s, &(m[i][j].v));
        }
    }
}

void set_cbcb(Vec *a, Vec *b, Tri **m, int l) {
    for (int i = 1; i <= l; i++) {
        Vec ai, bi, ci;
        vsub(a[i + 1], a[i - 1], &ai);
        vnorm(&ai);
        vsub(b[i], a[i], &bi);
        vnorm(&bi);
        vprod(ai, bi, &ci);
        
        for (int j = 1; j <= l; j++) {
            Vec aj, bj, cj;
            vsub(a[j + 1], a[j - 1], &aj);
            vnorm(&aj);
            vsub(b[j], a[j], &bj);
            vnorm(&bj);
            vprod(aj, bj, &cj);
            
            m[i][j].cos.x = vdot(ai, aj);
            m[i][j].cos.y = vdot(bi, bj);
            m[i][j].cos.z = vdot(ci, cj);
        }
    }
}

void setframe(Vec a, Vec b, Vec c, Mat *frame) {
    Vec x, y, z;
    
    vsub(c, a, &x);
    vave(c, a, &c);
    vsub(c, b, &y);
    vprod(y, x, &z);
    vprod(z, x, &y);
    vnorm(&x);
    vnorm(&y);
    vnorm(&z);
    VtoM(x, y, z, frame);
}

/* ============================================================================
 * MAIN FUNCTION
 * ============================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <pdb_file> [spread] [nruns] [subdom]\n", argv[0]);
        fprintf(stderr, "\nParameters:\n");
        fprintf(stderr, "  spread: Contact distance cutoff (default: 15)\n");
        fprintf(stderr, "  nruns:  Analysis mode (default: 1)\n");
        fprintf(stderr, "  subdom: Subdomain analysis flag (default: 0)\n");
        return 1;
    }
    
    DomainConfig config;
    if (!domain_analysis_init(&config, argc, argv)) {
        return 1;
    }
    
    srand(234);  /* Use standard srand instead of srand48 for portability */
    
    Pdbentry_ *prot = get_pdb(argv[1], 1, 1);
    if (!prot) {
        fprintf(stderr, "Error: Failed to load PDB file %s\n", argv[1]);
        return 1;
    }
    
    /* Check if we have valid data */
    if (prot->Chainno == 0) {
        fprintf(stderr, "Error: No chains found in PDB file\n");
        fprintf(stderr, "Note: This program only processes CA atoms. Check your PDB file contains CA atoms.\n");
        free_pdb(prot);
        return 1;
    }
    
    printf("Processing: %s\n", prot->Compound);
    printf("Found %d chain(s)\n", prot->Chainno);
    
    Seq seq = {0};
    Tri **triangles = NULL;
    int len = protin(prot, &seq, 0, &triangles, 1.0, 0);
    
    if (len == 0) {
        fprintf(stderr, "Error: Failed to process protein structure\n");
        free_pdb(prot);
        return 1;
    }
    
    printf("Processing %d residues\n", len);
    
    check_chain_breaks(&seq);
    
    if (config.subdom && len < 250) {
        printf("Only domains of 250 and over are considered for reparsing\n");
        free_pdb(prot);
        return 0;
    }
    
    int nn = len + 2;
    float **mat = alloc_matrix(nn);
    float **net = alloc_matrix(nn);
    float *dom = calloc(nn, sizeof(float));  /* Use calloc for initialization */
    int **doms = calloc(4, sizeof(int *));
    if (!doms) {
        fprintf(stderr, "Error: Memory allocation failed (doms rows)\n");
        free_matrix(mat, nn);
        free_matrix(net, nn);
        free(dom);
        if (triangles) {
            for (int i = 0; i <= len + 1; i++) free(triangles[i]);
            free(triangles);
        }
        free(seq.res); free(seq.acc); free(seq.dom); free(seq.rid); free(seq.ca); free(seq.cb);
        free_pdb(prot);
        return 1;
    }
    for (int k = 0; k < 4; k++) {
        doms[k] = calloc(len + 2, sizeof(int));
        if (!doms[k]) {
            fprintf(stderr, "Error: Memory allocation failed (doms[%d])\n", k);
            for (int t = 0; t < k; t++) free(doms[t]);
            free(doms);
            free_matrix(mat, nn);
            free_matrix(net, nn);
            free(dom);
            if (triangles) {
                for (int i = 0; i <= len + 1; i++) free(triangles[i]);
                free(triangles);
            }
            free(seq.res); free(seq.acc); free(seq.dom); free(seq.rid); free(seq.ca); free(seq.cb);
            free_pdb(prot);
            return 1;
        }
    }
    
    if (!mat || !net || !dom || !doms) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free_matrix(mat, nn);
        free_matrix(net, nn);
        free(dom);
        free_int_matrix(doms, 4);
        if (triangles) {
            for (int i = 0; i <= len + 1; i++) {
                free(triangles[i]);
            }
            free(triangles);
        }
        free(seq.res);
        free(seq.acc);
        free(seq.dom);
        free(seq.rid);
        free(seq.ca);
        free(seq.cb);
        free_pdb(prot);
        return 1;
    }
    
    int domn[4] = {0};
    
    
    if (config.nruns < 1) {
        for (int i = 1; i <= len; i++) {
            vcopy(seq.ca[i], seq.cb + i);
        }
        
        if (config.nruns == 0 || config.nruns == -2) {
            flatten(seq.cb, len, 5);
            config.spread += 3;
        }
        
        if (config.nruns > -2) {
            beta(net, mat, seq.ca, dom, len);
        } else {
            for (int i = 1; i <= len; i++) {
                dom[i] = (float)i;
            }
        }
        
        printf("Spread: %d\n", config.spread);
        for (int i = 1; i <= IT && define(0, i, &seq, seq.cb, dom, net, mat, len, domn, &config); i++);
        
        for (int i = 1; i <= len; i++) {
            doms[0][i] = seq.dom[i];
        }
        
        if (config.sheet) {
            float cuts = betacut(net, doms[0], len);
            printf("Error in beta split = %7.2f\n", cuts);
        }
    } else {
        printf("\n\n===============  normal  ===============\n\n");
        printf("Spread: %d\n", config.spread);
        
        if (config.nruns == 2) {
            beta(net, mat, seq.ca, dom, len);
        } else {
            for (int i = 1; i <= len; i++) {
                dom[i] = (float)i;
            }
        }
        
        for (int i = 1; i <= IT && define(0, i, &seq, seq.ca, dom, net, mat, len, domn + 0, &config); i++);
        for (int i = 1; i <= len; i++) {
            doms[0][i] = seq.dom[i];
        }
        
        printf("\n\n=============== smoothed ===============\n\n");
        for (int i = 1; i <= len; i++) {
            vcopy(seq.ca[i], seq.cb + i);
        }
        flatten(seq.cb, len, 5);
        config.spread += 3;
        printf("Spread: %d\n", config.spread);
        
        if (config.nruns == 2) {
            beta(net, mat, seq.ca, dom, len);
        } else {
            for (int i = 1; i <= len; i++) {
                dom[i] = (float)i;
            }
        }
        
        for (int i = 1; i <= IT && define(2, i, &seq, seq.cb, dom, net, mat, len, domn + 2, &config); i++);
        for (int i = 1; i <= len; i++) {
            doms[2][i] = seq.dom[i];
        }
        
        printf("\n\n");
        compare(net, doms[0], doms[2], domn[0], domn[2], len, config.sheet);
    }
    
    /* Cleanup */
    free_matrix(mat, nn);
    free_matrix(net, nn);
    free(dom);
    for (int k = 0; k < 4; k++) free(doms[k]);
    free(doms);
    
    if (triangles) {
        for (int i = 0; i <= len + 1; i++) {
            free(triangles[i]);
        }
        free(triangles);
    }
    
    free(seq.res);
    free(seq.acc);
    free(seq.dom);
    free(seq.rid);
    free(seq.ca);
    free(seq.cb);
    free_pdb(prot);
    
    printf("\nDomain analysis complete.\n");
    return 0;
}
