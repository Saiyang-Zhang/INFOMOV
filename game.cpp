#include "precomp.h"
#include "game.h"

#define LINES		1024
#define LINEFILE	"lines1024.dat"
#define ITERATIONS	16

#define IRand(x) ((int)(RandomFloat()*(x)))

int lx1[LINES], ly1[LINES], lx2[LINES], ly2[LINES];			// lines: start and end coordinates
int x1_, y1_, x2_, y2_;										// room for storing line backup
__int64 fitness = 0xfffffffff;								// similarity to reference image
int lidx = 0;												// current line to be mutated
float peak = 0;												// peak line rendering performance
Surface* reference, *backup;								// surfaces
int* ref8;													// grayscale image for evaluation
Timer tm;													// stopwatch

// -----------------------------------------------------------
// Mutate
// Randomly modify or replace one line.
// -----------------------------------------------------------
void MutateLine( int i )
{
	// backup the line before modifying it
	x1_ = lx1[i], y1_ = ly1[i];
	x2_ = lx2[i], y2_ = ly2[i];
	do
	{
		if (rand() & 1)
		{
			// small mutation (50% probability)
			lx1[i] += IRand( 6 ) - 3, ly1[i] += IRand( 6 ) - 3;
			lx2[i] += IRand( 6 ) - 3, ly2[i] += IRand( 6 ) - 3;
			// ensure the line stays on the screen
			lx1[i] = min( SCRWIDTH - 1, max( 0, lx1[i] ) );
			lx2[i] = min( SCRWIDTH - 1, max( 0, lx2[i] ) );
			ly1[i] = min( SCRHEIGHT - 1, max( 0, ly1[i] ) );
			ly2[i] = min( SCRHEIGHT - 1, max( 0, ly2[i] ) );
		}
		else
		{
			// new line (50% probability)
			lx1[i] = IRand( SCRWIDTH ), lx2[i] = IRand( SCRWIDTH );
			ly1[i] = IRand( SCRHEIGHT ), ly2[i] = IRand( SCRHEIGHT );
		}
	} while ((abs( lx1[i] - lx2[i] ) < 3) || (abs( ly1[i] - ly2[i] ) < 3));
}

void UndoMutation( int i )
{
	// restore line i to the backuped state
	lx1[i] = x1_, ly1[i] = y1_;
	lx2[i] = x2_, ly2[i] = y2_;
}

// -----------------------------------------------------------
// DrawWuLine
// Anti-aliased line rendering.
// Straight from: 
// https://www.codeproject.com/Articles/13360/Antialiasing-Wu-Algorithm
// -----------------------------------------------------------
void DrawWuLine( Surface *screen, int X0, int Y0, int X1, int Y1, uint clrLine )
{
    /* Make sure the line runs top to bottom */
    if (Y0 > Y1)
    {
        int Temp = Y0; Y0 = Y1; Y1 = Temp;
        Temp = X0; X0 = X1; X1 = Temp;
    }
    
    /* Draw the initial pixel, which is always exactly intersected by
    the line and so needs no weighting */
    screen->Plot( X0, Y0, clrLine );
    
    int XDir, DeltaX = X1 - X0;
    if( DeltaX >= 0 )
    {
        XDir = 1;
    }
    else
    {
        XDir   = -1;
        DeltaX = 0 - DeltaX; /* make DeltaX positive */
    }
    
    /* Special-case horizontal, vertical, and diagonal lines, which
    require no weighting because they go right through the center of
    every pixel */
    int DeltaY = Y1 - Y0;
    
    unsigned short ErrorAdj;
    unsigned short ErrorAccTemp, Weighting;
    
    /* Line is not horizontal, diagonal, or vertical */
    unsigned short ErrorAcc = 0;  /* initialize the line error accumulator to 0 */
    
    BYTE rl = GetRValue( clrLine );
    BYTE gl = GetGValue( clrLine );
    BYTE bl = GetBValue( clrLine );
    double grayl = rl * 0.299 + gl * 0.587 + bl * 0.114;
    
    /* Is this an X-major or Y-major line? */
    if (DeltaY > DeltaX)
    {
    /* Y-major line; calculate 16-bit fixed-point fractional part of a
    pixel that X advances each time Y advances 1 pixel, truncating the
        result so that we won't overrun the endpoint along the X axis */
        ErrorAdj = ((unsigned long) DeltaX << 16) / (unsigned long) DeltaY;
        /* Draw all pixels other than the first and last */
        while (--DeltaY) {
            ErrorAccTemp = ErrorAcc;   /* remember currrent accumulated error */
            ErrorAcc += ErrorAdj;      /* calculate error for next pixel */
            if (ErrorAcc <= ErrorAccTemp) {
                /* The error accumulator turned over, so advance the X coord */
                X0 += XDir;
            }
            Y0++; /* Y-major, so always advance Y */
                  /* The IntensityBits most significant bits of ErrorAcc give us the
                  intensity weighting for this pixel, and the complement of the
            weighting for the paired pixel */
            Weighting = ErrorAcc >> 8;
            
            COLORREF clrBackGround = screen->pixels[X0 + Y0 * SCRWIDTH];
            BYTE rb = GetRValue( clrBackGround );
            BYTE gb = GetGValue( clrBackGround );
            BYTE bb = GetBValue( clrBackGround );
            double grayb = rb * 0.299 + gb * 0.587 + bb * 0.114;
            
            BYTE rr = ( rb > rl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( rb - rl ) + rl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( rl - rb ) + rb ) ) );
            BYTE gr = ( gb > gl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( gb - gl ) + gl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( gl - gb ) + gb ) ) );
            BYTE br = ( bb > bl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( bb - bl ) + bl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( bl - bb ) + bb ) ) );
            screen->Plot( X0, Y0, RGB( rr, gr, br ) );
            
            clrBackGround = screen->pixels[X0 + XDir + Y0 * SCRWIDTH];
            rb = GetRValue( clrBackGround );
            gb = GetGValue( clrBackGround );
            bb = GetBValue( clrBackGround );
            grayb = rb * 0.299 + gb * 0.587 + bb * 0.114;
            
            rr = ( rb > rl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( rb - rl ) + rl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( rl - rb ) + rb ) ) );
            gr = ( gb > gl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( gb - gl ) + gl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( gl - gb ) + gb ) ) );
            br = ( bb > bl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( bb - bl ) + bl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( bl - bb ) + bb ) ) );
            screen->Plot( X0 + XDir, Y0, RGB( rr, gr, br ) );
        }
        /* Draw the final pixel, which is always exactly intersected by the line
        and so needs no weighting */
        screen->Plot( X1, Y1, clrLine );
        return;
    }
    /* It's an X-major line; calculate 16-bit fixed-point fractional part of a
    pixel that Y advances each time X advances 1 pixel, truncating the
    result to avoid overrunning the endpoint along the X axis */
    ErrorAdj = ((unsigned long) DeltaY << 16) / (unsigned long) DeltaX;
    /* Draw all pixels other than the first and last */
    while (--DeltaX) {
        ErrorAccTemp = ErrorAcc;   /* remember currrent accumulated error */
        ErrorAcc += ErrorAdj;      /* calculate error for next pixel */
        if (ErrorAcc <= ErrorAccTemp) {
            /* The error accumulator turned over, so advance the Y coord */
            Y0++;
        }
        X0 += XDir; /* X-major, so always advance X */
                    /* The IntensityBits most significant bits of ErrorAcc give us the
                    intensity weighting for this pixel, and the complement of the
        weighting for the paired pixel */
        Weighting = ErrorAcc >> 8;
        
        COLORREF clrBackGround = screen->pixels[X0 + Y0 * SCRWIDTH];
        BYTE rb = GetRValue( clrBackGround );
        BYTE gb = GetGValue( clrBackGround );
        BYTE bb = GetBValue( clrBackGround );
        double grayb = rb * 0.299 + gb * 0.587 + bb * 0.114;
        
        BYTE rr = ( rb > rl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( rb - rl ) + rl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( rl - rb ) + rb ) ) );
        BYTE gr = ( gb > gl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( gb - gl ) + gl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( gl - gb ) + gb ) ) );
        BYTE br = ( bb > bl ? ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( bb - bl ) + bl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?Weighting:(Weighting ^ 255)) ) / 255.0 * ( bl - bb ) + bb ) ) );
        
        screen->Plot( X0, Y0, RGB( rr, gr, br ) );
        
        clrBackGround = screen->pixels[X0 + (Y0 + 1 )* SCRWIDTH];
        rb = GetRValue( clrBackGround );
        gb = GetGValue( clrBackGround );
        bb = GetBValue( clrBackGround );
        grayb = rb * 0.299 + gb * 0.587 + bb * 0.114;
        
        rr = ( rb > rl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( rb - rl ) + rl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( rl - rb ) + rb ) ) );
        gr = ( gb > gl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( gb - gl ) + gl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( gl - gb ) + gb ) ) );
        br = ( bb > bl ? ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( bb - bl ) + bl ) ) : ( ( BYTE )( ( ( double )( grayl<grayb?(Weighting ^ 255):Weighting) ) / 255.0 * ( bl - bb ) + bb ) ) );
        
        screen->Plot( X0, Y0 + 1, RGB( rr, gr, br ) );
    }
    
    /* Draw the final pixel, which is always exactly intersected by the line
    and so needs no weighting */
    screen->Plot( X1, Y1, clrLine );
}

// -----------------------------------------------------------
// Fitness evaluation
// Compare current generation against reference image.
// -----------------------------------------------------------
__int64 Game::Evaluate()
{
	// compare to reference using SIMD magic. don't worry about it, it's fast.
	uint* src = screen->pixels;
	const int quads = (SCRWIDTH * SCRHEIGHT) / 4;
	__m128i* A4 = (__m128i*)src;
	__m128i* B4 = (__m128i*)ref8;
	union { __m128i diff4; int diff[4]; };
	diff4 = _mm_set1_epi32( 0 );
	union { __m128i mask4; int mask[4]; };
	mask[0] = mask[1] = mask[2] = mask[3] = 255;
	for (int i = 0; i < quads; i++)
	{
		const __m128i d2 = _mm_abs_epi32( _mm_sub_epi32( _mm_and_si128( A4[i], mask4 ), B4[i] ) );
		diff4 = _mm_add_epi32( diff4, _mm_srai_epi32( _mm_mul_epi32( d2, d2 ), 12 ) );
	}
	__int64 retval = diff[0];
	retval += diff[1];
	retval += diff[2];
	retval += diff[3];
	return retval;
}

// -----------------------------------------------------------
// Application initialization
// Load a previously saved generation, if available.
// -----------------------------------------------------------
void Game::Init()
{
	for (int i = 0; i < LINES; i++) MutateLine( i );
	FILE* f = fopen( LINEFILE, "rb" );
	if (f)
	{
		fread( lx1, 4, LINES, f );
		fread( ly1, 4, LINES, f );
		fread( lx2, 4, LINES, f );
		fread( ly2, 4, LINES, f );
		fclose( f );
	}
	Surface* reference = new Surface( "assets/image3.png" );
	backup = new Surface( SCRWIDTH, SCRHEIGHT );
	ref8 = (int*)MALLOC64( SCRWIDTH * SCRHEIGHT * 4 );
	for (int i = 0; i < (SCRWIDTH * SCRHEIGHT); i++) ref8[i] = reference->pixels[i] & 255;
	fitness = 512 * 512 * 16;
}

// -----------------------------------------------------------
// Application termination
// Save the current generation, so we can continue later.
// -----------------------------------------------------------
void Game::Shutdown()
{
	FILE* f = fopen( LINEFILE, "wb" );
	fwrite( lx1, 4, LINES, f );
	fwrite( ly1, 4, LINES, f );
	fwrite( lx2, 4, LINES, f );
	fwrite( ly2, 4, LINES, f );
	fclose( f );
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float _DT )
{
	tm.reset();
	int lineCount = 0;
	int iterCount = 0;
	// draw up to lidx
	memset( screen->pixels, 255, SCRWIDTH * SCRHEIGHT * 4 );
	for (int j = 0; j < lidx; j++, lineCount++)
	{
		unsigned int c = (j * 128) / LINES;
		DrawWuLine( screen, lx1[j], ly1[j], lx2[j], ly2[j], c + (c << 8) + (c << 16) );
	}
	int base = lidx;
	screen->CopyTo( backup, 0, 0 );
	// iterate and draw from lidx to end
	for (int k = 0; k < ITERATIONS; k++)
	{
		memcpy( screen->pixels, backup->pixels, SCRWIDTH * SCRHEIGHT * 4 );
		MutateLine( lidx );
		for (int j = base; j < LINES; j++, lineCount++)
		{
			unsigned int c = (j * 128) / LINES;
			DrawWuLine( screen, lx1[j], ly1[j], lx2[j], ly2[j], c + (c << 8) + (c << 16) );
		}
		__int64 diff = Evaluate();
		if (diff < fitness) fitness = diff; else 
		UndoMutation( lidx );
		// lidx = (lidx + 1) % LINES;
		iterCount++;
	}
	// stats
	char t[128];
	float elapsed = tm.elapsed();
	float lps = (float)lineCount / elapsed;
	peak = max( lps, peak );
	sprintf( t, "fitness: %i", fitness );
	screen->Bar( 0, SCRHEIGHT - 33, 130, SCRHEIGHT - 1, 0 );
	screen->Print( t, 2, SCRHEIGHT - 24, 0xffffff );
	sprintf( t, "lps:     %5.2fK", lps );
	screen->Print( t, 2, SCRHEIGHT - 16, 0xffffff );
	sprintf( t, "ips:     %5.2f", (iterCount * 1000) / elapsed );
	screen->Print( t, 2, SCRHEIGHT - 8, 0xffffff );
	sprintf( t, "peak:    %5.2f", peak );
	screen->Print( t, 2, SCRHEIGHT - 32, 0xffffff );
}