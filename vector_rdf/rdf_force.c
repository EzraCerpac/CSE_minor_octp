
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINE 1024
#define BinNum 3000
#define Delta 0.005

int parse_dump_timestep(FILE *stream, long *Timestep, long *NumberOfParticles, 
                        double *Box, long n, double Rxx[], double Ryy[], 
                        double Rzz[], double Fxx[], double Fyy[], double Fzz[]) 
{
    /**
     * @brief Parse a single timestep from a LAMMPS dump file created with 
                    "custom ${Nrun} dump id x y z fx fy fz"
     * @param stream File pointer of the dump file with its offset at the 
                        beginning of the timestep, offset is always moved to 
                        the beginning of the next timestep
     * @param NumberOfParticles Return addres for the number of particles in the timestep
     * @param Box Return address for the box size, box is assumed to be cubic 
                    otherwise the function exits and returns 1
     * @param n Size of the arrays Rxx, Ryy, Rzz, Fxx, Fyy, Fzz, 
                    if the number of particles in the timestep is larger than 
                    n, the function exits and returns 2
     * @param Rxx Return array for the x coordinates of the particles
     * @param Ryy Return array for the y coordinates of the particles
     * @param Rzz Return array for the z coordinates of the particles
     * @param Fxx Return array for the x components of the forces on the particles
     * @param Fyy Return array for the y components of the forces on the particles
     * @param Fzz Return array for the z components of the forces on the particles
     * @return 0 if the timestep was parsed successfully, 1 if an error occured
     */
    char buffer[MAXLINE];
    double boxLeft[3], boxRight[3];
    long i, ID;

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: TIMESTEP") != NULL) {
        fgets(buffer, MAXLINE, stream);
        *Timestep = atoi(buffer);
    } else if (feof(stream)) {
        return 1;
    } else {
    fprintf(stderr, "Expected next line to be 'ITEM: TIMESTEP', but instead \
                    is:\n%s.\n", buffer);
    return 1;   
    }

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: NUMBER OF ATOMS") != NULL) {
        fgets(buffer, MAXLINE, stream);
        *NumberOfParticles = atoi(buffer);
        if (*NumberOfParticles > n) {
            fprintf(stderr, "Size of the return arrays are not large enough to \
            hold the number of particles : %ld.\n", *NumberOfParticles);
            return 1;
        }
    } else {
    fprintf(stderr, "Expected next line to be 'ITEM: NUMBER OF ATOMS', but \
                    instead is:\n%s.\n", buffer);
    return 1;
    }

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: BOX BOUNDS") != NULL) {
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[0], &boxRight[0]);
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[1], &boxRight[1]);
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[2], &boxRight[2]);
        // TODO: this is comparing floats, maybe change to checking "if close"
        if (boxLeft[0] != boxLeft[1] || boxLeft[0] != boxLeft[2] || 
            boxRight[0] != boxRight[1] || boxRight[0] != boxRight[2]) {
            fprintf(stderr, "Box is not a cube.\n");
            return 1;
        }
    } else {
    fprintf(stderr, "Expected next line to be 'ITEM: BOX BOUNDS', \
                    but instead is:\n%s.\n", buffer);
    return 1;
    }

    *Box = boxRight[0] - boxLeft[0];

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: ATOMS") == NULL) {
    fprintf(stderr, "Expected next line to be 'ITEM: ATOMS', but instead \
                    is:\n%s.\n", buffer);
    return 1;
    } else{
        for (i = 0; i < *NumberOfParticles; i++) {
            fgets(buffer, MAXLINE, stream);
            sscanf(buffer, "%ld %lf %lf %lf %lf %lf %lf", &ID, &Rxx[i], 
                    &Ryy[i], &Rzz[i], &Fxx[i], &Fyy[i], &Fzz[i]);
        }
    }
    return 0;
}


void sample_rdf(long NumberOfParticles, double Box, double g[], double Rxx[], 
                double Ryy[], double Rzz[], double Fxx[], double Fyy[], 
                double Fzz[]) 
{
    /**
     * @brief Sample the radial distribution function of the system in array g 

     * @param NumberOfParticles Number of particles in the system
     * @param Box Size of the box
     * @param g Array to store the sampled rdf, array values are zeroed before sampling.
     * @param Rxx Array of x coordinates of the particles
     * @param Ryy Array of y coordinates of the particles
     * @param Rzz Array of z coordinates of the particles
     * @param Fxx Array of x components of the forces on the particles
     * @param Fyy Array of y components of the forces on the particles
     * @param Fzz Array of z components of the forces on the particles
     */
    long i, j, k;
    double Dx, Dy, Dz, Rij, Ff, R;
    double HalfBox = Box/2.0;

    for (i = 0; i < BinNum; i++) {
        g[i] = 0.0;
    }

    for(i=0;i<NumberOfParticles;i++) {
        for(j=i+1;j<NumberOfParticles;j++) {
            
            // Calculate the distance between the two particles, R, accounting
            // for periodic boundary conditions
            Dx = Rxx[i] - Rxx[j];
            Dy = Ryy[i] - Ryy[j];
            Dz = Rzz[i] - Rzz[j];
 
            if (Dx > HalfBox) {
               Dx = Dx - Box;
            } else if (Dx < -HalfBox) {
               Dx = Dx + Box;
            }
 
            if (Dy > HalfBox) {
               Dy = Dy - Box;
            } else if (Dy< -HalfBox) {
               Dy = Dy + Box;
            }
 
            if (Dz > HalfBox) {
               Dz = Dz - Box;
            } else if (Dz < -HalfBox) {
               Dz = Dz + Box;
            }

            Rij = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

            // RDF is only valid till half the box size
            if (Rij < HalfBox) {
                Ff = (Fxx[i]*Dx + Fyy[i]*Dy + Fzz[i]*Dz) / (Rij*Rij*Rij);

                for (k = 0; k < BinNum; k++) {
                    R = (k + 0.5) * Delta;

                if (R < HalfBox && R > Rij)
                    g[k] += Ff;
                }
            }
        }
    }
}

int main()
{
    FILE *fp;
    int status;
    long NumberOfParticles, i, k, n, timestep, nTimesteps, minTimestep;
    double Temp, c, Box, R, *Rxx, *Ryy, *Rzz, *Fxx, *Fyy, *Fzz, *g, *gavg;

    minTimestep = 20000;
    Temp = 299;
    n = 1728;
    
    fp = fopen("dump", "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open file 'dump'\n");
        return 1;
    }

    Rxx = (double *) calloc(n, sizeof(double));
    Ryy = (double *) calloc(n, sizeof(double));
    Rzz = (double *) calloc(n, sizeof(double));
    Fxx = (double *) calloc(n, sizeof(double));
    Fyy = (double *) calloc(n, sizeof(double));
    Fzz = (double *) calloc(n, sizeof(double));
    g = (double *) malloc(BinNum * sizeof(double));
    gavg = (double *) calloc(BinNum, sizeof(double));    

    nTimesteps = 0;
    while (feof(fp) == 0) {
        status = parse_dump_timestep(fp, &timestep, &NumberOfParticles, &Box, n, 
                                    Rxx, Ryy, Rzz, Fxx, Fyy, Fzz);

        if ((status == 0) && (timestep >= minTimestep)) {
            sample_rdf(NumberOfParticles, Box, g, Rxx, Ryy, Rzz, Fxx, Fyy, Fzz);
            c = Box*Box*Box / (2 * M_PI * Temp * NumberOfParticles*NumberOfParticles); 
            // TODO: read the correct Temp for each timestep from temperature.dat
            for (k = 0; k < BinNum; k++) {
                gavg[k] += g[k] * c;
            }
            nTimesteps++;
        }
    }

    if (nTimesteps != 0) {
        for (k = 0; k < BinNum; k++) {
            gavg[k] /= nTimesteps; // Average the rdf
        }
    }
    printf("nTimesteps = %ld\n", nTimesteps);
    fclose(fp);

    // Write the rdf to a file
    fp = fopen("frdf.dat", "w");
    for (i = 0; i < BinNum; i++) {
        R = (i + 0.5) * Delta;
        if (R < Box/2.0){
            fprintf(fp, "%lf %lf\n", R, gavg[i]);
        }
    }
    fclose(fp);
    return 0;
}