/* sbox.c   
 *
 * by: David Canright
 *
 * illustrates compact implementation of AES S-box via subfield operations
 *   case # 4 : [d^16, d], [alpha^8, alpha^2], [Omega^2, Omega]
 *   nu = beta^8 = N^2*alpha^2, N = w^2
 */

#include <stdio.h>
#include <sys/types.h>
#include "./matrix.h"

//static int Clefia_transofrm[8]= {0x0A,0x41,0x58,0x20,0x30,0x02,0x90,0x44};
//static int Clefia_transofrm[8]= {0x44,0x90,0x02,0x30,0x20,0x58,0x41,0x0A};
//static int Clefia_transofrm[8]= {0x1e,0x51,0x01,0x06,0x65,0x5A,0x60,0x81};
//static int Clefia_transofrm[8]= {0x81,0x60,0x5A,0x65,0x06,0x01,0x51,0x1e};
//static int Clefia_transofrm[8]= {0x18,0x51,0x01,0x06,0x65,0x5c,0x60,0x81};
//static int Clefia_transofrm[8]= {0x18,0x88,0x80,0x60,0xA6,0x3A,0x06,0x81};
static int Clefia_transofrm[8]= {0x01,0x4e,0x0A,0xc4,0x84,0x1c,0x10,0x69};

#define CIPHER 0

#define CONSTANT 0x69
//#define CONSTANT 0x63

/* multiply in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_mul( int x, int y ) {
    int a, b, c, d, e, p, q;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    c = (y & 0x2) >> 1; d = (y & 0x1);
    e = (a ^ b) & (c ^ d);
    p = (a & c) ^ e;
    q = (b & d) ^ e;
    return ( (p<<1) | q );
    }

/* scale by N = Omega^2 in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N( int x ) { //Figure 16
    int a, b, p, q;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    p = b;
    q = a ^ b;
    return ( (p<<1) | q );
    }

/* scale by N^2 = Omega in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N2( int x ) { //Figure 15
    int a, b, p, q;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    p = a ^ b;
    q = a;
    return ( (p<<1) | q );
    }

/* scale by N^2 = Omega^2 in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N2_w2_scaler( int x ) { //Figure 18
    int a, b, p, q;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    p = a;
    q = a^ b;
    return ( (p<<1) | q );
    }
/* scale by N = Omega in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N_w_scaler( int x ) {  //Figure 17
    int a, b, p, q;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    p = a ^ b;
    q = b;
    return ( (p<<1) | q );
    }
/* square in GF(2^2), using normal basis (Omega^2,Omega) */
/* NOTE: inverse is identical */
int G4_sq( int x ) {
    int a, b;

    a = (x & 0x2) >> 1; b = (x & 0x1);
    return ( (b<<1) | a );
    }

/* multiply in GF(2^4), using normal basis (alpha^8,alpha^2) */
int G16_mul( int x, int y ) {
    int a, b, c, d, e, p, q;

    a = (x & 0xC) >> 2; b = (x & 0x3);
    c = (y & 0xC) >> 2; d = (y & 0x3);
    e = G4_mul( a ^ b, c ^ d );
    e = G4_scl_N(e);
    p = G4_mul( a, c ) ^ e;
    q = G4_mul( b, d ) ^ e;
    return ( (p<<2) | q );
    }

/* square & scale by nu in GF(2^4)/GF(2^2), normal basis (alpha^8,alpha^2) */
/*   nu = beta^8 = N^2*alpha^2, N = w^2 */
int G16_sq_scl( int x ) {
    int a, b, p, q;

    a = (x & 0xC) >> 2; b = (x & 0x3);
    // d16,d
    if(OPTION==0){
        // C=N, D=0
        p=G4_scl_N(G4_sq(a));
        q = G4_scl_N2(G4_sq(a^b));
    }
    else if(OPTION==1){

        // C=0, D=N
        p= G4_scl_N2(G4_sq(a^b)); 
        q= G4_scl_N(G4_sq(b));
    }

    else if(OPTION==2){
        // C=N2, D=0
        p=G4_scl_N2(G4_sq(a));
        q =G4_sq(a^b);
    }
    else if(OPTION==3){
        // C=0, D=N2
        p= (G4_sq(a^b)); 
        q= G4_scl_N2(G4_sq(b));
    }
    else if(OPTION==4){
        // C=N, D=1
        p= G4_scl_N(G4_sq(b)); 
        q= G4_scl_N2(G4_sq(a))^p;
    }
    else if(OPTION==5){
        // C=1, D=N
        q= G4_scl_N(G4_sq(a)); 
        q= G4_scl_N2(G4_sq(a))^q;
    }
    else if (OPTION==6){
        // C=N2, D=1
        p=G4_sq(a)^ G4_scl_N(G4_sq(b));
        q = G4_sq(a);
    }
    else if (OPTION==7){
        // C=1, D=N2
        p= G4_sq(b);
        q = G4_sq(b)^ G4_scl_N(G4_sq(a)); 
    }

    else{
    // d64,d4
     perror("not in the option list");
    }
    
    return ( (p<<2) | q );
    }

/* inverse in GF(2^4), using normal basis (alpha^8,alpha^2) */
int G16_inv( int x ) {
    int a, b, c, d, e, p, q;

    a = (x & 0xC) >> 2; b = (x & 0x3);
    c = G4_scl_N( G4_sq( a ^ b ) );
    d = G4_mul( a, b );
    e = G4_sq( c ^ d );   // really inverse, but same as square
    p = G4_mul( e, b );
    q = G4_mul( e, a );
    return ( (p<<2) | q );
    }

/* inverse in GF(2^8), using normal basis (d^16,d) */
int G256_inv( int x ) {
    int a, b, c, d, e, p, q;

    a = (x & 0xF0) >> 4; b = (x & 0x0F);
    c = G16_sq_scl( a ^ b );
    d = G16_mul( a, b );
    e = G16_inv( c ^ d );
    p = G16_mul( e, b );
    q = G16_mul( e, a );
    return ( (p<<4) | q );
    }

/* convert to new basis in GF(2^8) */
/* i.e., bit matrix multiply */
int G256_newbasis( int x, int b[] ) {
    int i, y = 0;

    for ( i=7; i >= 0; i-- ) {
        if ( x & 1 ) y ^= b[i];
        x >>= 1;
        }
    return ( y );
    }

/* find Sbox of n in GF(2^8) mod POLY */
int Sbox( int n ) {
     int t;

    t = G256_newbasis( n, A2X );
    t = G256_inv( t );
    t = G256_newbasis( t, X2S );
    return ( t ^ CONSTANT );
    }

/* find inverse Sbox of n in GF(2^8) mod POLY */
int iSbox( int n ) {
     int t;

    t = G256_newbasis( n ^ CONSTANT, S2X );
    t = G256_inv( t );
    t = G256_newbasis( t, X2A );
    return ( t );
    }



/* compute tables of Sbox & its inverse; print 'em out */
int main () {
     int Sbox_tbl[256], iSbox_tbl[256], i, j;
    int T;
    for (i = 0; i < 256; i++) {
        if (CIPHER == 0){ 
        T=G256_newbasis(i,Clefia_transofrm)^0x1e;    
        Sbox_tbl[i] = Sbox(T);
        //Sbox_tbl[i]=G256_newbasis(i,Clefia_transofrm)^0x1e;
        //Sbox_tbl[i] = G256_newbasis(Sbox(i),Clefia_transofrm)^0x69;       
        //printf("%02x\n",T);      
        iSbox_tbl[i] = 0;
        }

        else if( CIPHER == 1){
        Sbox_tbl[i] = Sbox(i);
        iSbox_tbl[i] = iSbox(i);
        }

        }
    printf ("char S[256] = {\n");
    for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
        printf ( "%02x, ", Sbox_tbl[i*16+j]);
        }
        printf ( "\n" );
        }
        printf ( "};\n\n" );
    printf ("char Si[256] = {\n");
    for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
        printf ( "%02x, ", iSbox_tbl[i*16+j]);
        }
        printf ( "\n" );
        }
        printf ( "};\n\n" );
    return(0);
    }