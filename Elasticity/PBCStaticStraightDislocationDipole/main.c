#include <NIST.h>
/* We are going to adopt the eV/A system */
#include <Atoms.h>
#include <Min.h>
#include <VecMat2.h>
#include <Elasticity.h>

Aapp_Define_Config;
Chemtab ct[1]={{0}};
Tp *tp=NULL;
Neighborlist N[1]={{0}};
bool save_screen_output;
double separation_removal;
int option;
ElasticityPBCStaticStraightDislocationDipole dipole[1];
Elasticity e = UnknownElasticity;
double dipole_s0[2], dipole_s1[2], dipole_b[4];
ElasticityPBCStaticStraightDislocationDipoleUserSetup user[1]=
{UnknownElasticityPBCStaticStraightDislocationDipoleUserSetup};

TermString fname_dipole;
TermString fname_tee;
TermString fname_perfect;

#define ANGULAR_MESH 180
TermString fname_angular="angular.out";

int main (int argc, char *argv[])
{
    int i;
    TermString buf;
    double *c, d;
    FILE *fp=NULL, *fp_angular;

    fscanf (stdin, "Output CFG file (default=dipole):\n%s\n\n",
            fname_dipole);
    if (!strcasecmp(fname_dipole, "default"))
        strcpy(fname_dipole, "dipole");

    fscanf (stdin, "Save screen output to file (default=*.out, "
            "NULL=don't):\n%s\n\n", fname_tee);
    if ( (save_screen_output=strcasecmp(fname_tee, "NULL")) )
    {
        if ( !strcasecmp(fname_tee, "default") )
            sprintf (fname_tee, "%s.out", fname_dipole);
        fp = wopen (fname_tee);
        ftie (stdout, fp);
    }
    else ft = stdout;
    
    fscanf (stdin, "CFG file of a stress-free crystal "
            "(default=perfect_crystal):\n%s\n\n", fname_perfect);
    if (!strcasecmp(fname_perfect, "default"))
        sprintf(fname_perfect, "%s", "perfect_crystal");
    Config_Load (fname_perfect, ft, Config_Aapp_to_Alib);
    Config_TRIM(Config_Aapp_to_Alib);

    fscanf (stdin, "Elastic constant matrix in the frame CFG is "
            "presented [GPa]:\n"
            "C11     C12     C13     C14     C15     C16\n"
            "C21     C22     C23     C24     C25     C26\n"
            "C31     C32     C33     C34     C35     C36\n"
            "C41     C42     C43     C44     C45     C46\n"
            "C51     C52     C53     C54     C55     C56\n"
            "C61     C62     C63     C64     C65     C66\n"
            "---------------------------------------------"
            "------------------------\n"
            "%lf %lf %lf %lf %lf %lf\n" 
            "%lf %lf %lf %lf %lf %lf\n" 
            "%lf %lf %lf %lf %lf %lf\n" 
            "%lf %lf %lf %lf %lf %lf\n" 
            "%lf %lf %lf %lf %lf %lf\n" 
            "%lf %lf %lf %lf %lf %lf\n\n",
            M6e(e.C));
    M6DividE(e.C, USTRESS_IN_GPA);
    S6fPRmul (ft, "C = %M GPa\n ", e.C, USTRESS_IN_GPA);

    fscanf (stdin, "Burgers vector of the first dislocation "
            "(whose xi is // the third edge\n"
            "of the PBC box) in the frame CFG is presented [A]:\n"
            "%lf %lf %lf\n\n", V3e(dipole_b));
    dipole_b[3] = V3LENGTH(dipole_b);
    V3fpr(ft,    "Burgers vector = %M A\n", dipole_b);
    fprintf (ft, "        length = %g A\n\n", dipole_b[3]);

    fscanf (stdin, "First dislocation reduced coordinates:\n"
            "%lf %lf\n\n", V2e(dipole_s0));
    V2fpr(ft, "dipole_s0 = %M\n", dipole_s0);

    fscanf (stdin, "Second dislocation reduced coordinates:\n"
            "%lf %lf\n\n", V2e(dipole_s1));
    V2fpr(ft, "dipole_s1 = %M\n ", dipole_s1);

    fscanf (stdin, "Remove atoms closer than [A] "
            "(default=0.5, negative means don't):\n%s\n\n", buf);
    if (!strcasecmp(buf, "default")) separation_removal = 0.5;
    else sscanf(buf, "%lf", &separation_removal);
    fprintf (ft, "Atoms closer than %g A will be removed;\n\n",
             separation_removal);

    fscanf (stdin, "Control option (choose 0 to 4, default=4):\n%s\n\n", buf);
    if (!strcasecmp(buf, "default")) option = 4;
    else sscanf(buf, "%d", &option);
    switch(option)
    {
        case 0: c = &(user->H[0][0]); break;
        case 1: c = &(user->D[0][0]); break;
        case 2: c = &(user->Stress0[0][0]); break;
        case 3: c = &(user->Stress1[0][0]); break;
        case 4: c = &(user->Stresscell[0][0]); break;
        default:
            pe("option = %d outside of 0..4 range\n", option);
    }
    for (i=0; i<=4; i++)
    {
        switch(i)
        {
            case 0: fscanf (stdin, "(Option 0) "
                            "H-matrix (new PBC box) [A]:\n");
            break;
            case 1: fscanf (stdin, "(Option 1) "
                            "D-matrix (H=H0*(1+D), strain=(D+D^T)/2):\n");
            break;
            case 2: fscanf (stdin, "(Option 2) "
                            "Local stress at the first dislocation [GPa]:\n");
            break;
            case 3: fscanf (stdin, "(Option 3) "
                            "Local stress at the second dislocation [GPa]:\n");
            break;
            case 4: fscanf (stdin, "(Option 4) "
                            "Supercell Virial stress average [GPa]:\n");
            break;
        }
        if (i == option)
        {
            fscanf (stdin, "%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n\n",
                    V9e(c));
            if ( (i==2) || (i==3) || (i==4) )
                V9DiV(c, USTRESS_IN_GPA);
        }
        else fscanf (stdin, "%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n\n",
                     &d, &d, &d, &d, &d, &d, &d, &d, &d);
    }

    ElasticityPBCStaticStraightDislocationDipoleASSIGN
        (H, e.C, dipole_s0, dipole_s1, dipole_b, dipole);

    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, user);

    S3fPR(ft, "H = %M A\n ", user->H);
    S3fPR(ft, "D = %M A\n ", user->D);
    S3fPR(ft, "simpleStrain = %M A\n ", user->simpleStrain);
    S3fPRmul(ft, "Stress0 = %M GPa\n ", user->Stress0, USTRESS_IN_GPA);
    V3fpr(ft, "Force0 = %M eV/A\n ", user->Force0);
    S3fPRmul(ft, "Stress1 = %M GPa\n ", user->Stress1, USTRESS_IN_GPA);
    V3fpr(ft, "Force1 = %M eV/A\n ", user->Force1);
    S3fPRmul(ft, "Stresscell = %M GPa\n ", user->Stresscell, USTRESS_IN_GPA);

    fprintf (ft, "Ebench = %g eV.\n\n"
             "Remember:\n"
             "Eatomistic = Ebench + 2*(Ecore-PE*log(r0[A]))*linelength;\n"
             "here, PE = %g eV/A, linelength = %g A.\n\n",
             user->Ebench, dipole->d->PE, dipole->linelength);

    for (i=np; i--;)
        ElasticityPBCStaticStraightDislocationDipoleTRANSFORM
            (dipole, s+3*i);
    M3EQV(user->H, H);

    Config_SAVE (Config_Aapp_to_Alib, fname_dipole);
    fprintf (ft, "The dipole configuration has been saved on \"%s\".\n\n",
             fname_dipole);

    fprintf (ft, "Note: the core energy definition depends "
             "on the definition of theta = 0,\n");
    fprintf (ft, "which is taken to be (%.5f, %.5f, %.5f) direction now,\n"
             "for which A = 0 [eV/A].\n", V3E(dipole->d->m));
    fprintf (ft, "For non-zero theta, the core energy needs to be added "
             "by A(theta) [eV/A].\n");
    fp_angular = wopen (fname_angular);
    for (i=0; i<=ANGULAR_MESH; i++)
    {
        d = PI / ANGULAR_MESH * i;
        fprintf (fp_angular, "%g %g\n", RADIAN_TO_DEGREE(d),
                 ElasticityPBCStaticStraightDislocationDipoleAngularEnergy
                 (dipole, d));
    }
    fclose (fp_angular);
    fprintf (ft, "A(theta) [eV/A] vs theta [degree] has been "
             "saved on \"%s\".\n\n", fname_angular);
    
    fbrk();
    fclose(fp);
    if (save_screen_output)
        printf("The screen output has been saved on \"%s\".\n\n", fname_tee);

    printf("Mission succeed.\n\n");
    return(0);
} /* end main() */
