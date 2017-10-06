/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/* Compute CNA and assign it as auxilary   */
/* Thu Apr 11 2013 Ju Li <liju@mit.edu>    */
/*******************************************/


#include "A.h"

int CNA_hist[5];
const Atom_CNA_color ATOM_CNA_COLOR
[5] =
{
    {"Grey",             0.6,      0.6,      0.6       },
    {"Red",               1.,       0.,        0.      },
    {"Green",             0.,       1.,        0.      },
    {"Blue",              0.,       0.,        1.      },
    {"Yellow",          0.95,0.7900735, 0.0138587      },
};


int ComputeCNA = 0;
short *CNA_pattern=NULL;

/* Free the auxiliary properties based on CNA */
void CNA_Free()
{
    ComputeCNA = 0;
    Free(CNA_pattern);
    return;
} /* end CNA_Free() */

/* Append CNA as auxiliary properties */
void CNA_ListBuild()
{
    
    int i,iatom;
    short pattern=0;
    for (i=0; i<5; i++) CNA_hist[i] = 0;
    
    Neighborlist N_NonPairwise[1];
    Neighborlist_create_nonpairwise_compressed_image
    (Config_Aapp_to_Alib, ct,&tp, N, N_NonPairwise);
    
    REALLOC (CNA_ListBuild, CNA_pattern, np, short);
    
    for (iatom=0; iatom<np; iatom++){
        
        int j,k,l,jatom,katom,NoCommon,NoNeighs,maxbond,minbond,nobonds;
        int CommonNeighList[14],CommonBonds[14];
        NoNeighs = coordination[iatom];
        
        if (NoNeighs!=12 && NoNeighs!=14) pattern=CNA_OTHER;
        else if (NoNeighs==12){
            int nfcc,nhcp,nico;
            nfcc=nhcp=nico=0;
            for (j=N_NonPairwise->idx[iatom];
                 j<N_NonPairwise->idx[iatom+1]; j++){
                NoCommon = 0;
                maxbond = 0;
                minbond = ATOM_COORDINATION_MAX;
                nobonds = 0;
                jatom = N_NonPairwise->list[j];
                
                for (k=N_NonPairwise->idx[iatom];
                     k<N_NonPairwise->idx[iatom+1]; k++){
                    katom = N_NonPairwise->list[k];
                    if(CNA_CheckBond(jatom,katom))
                        CommonNeighList[NoCommon++] = katom;
                }
                
                if (NoCommon==4 || NoCommon==5){
                    for (k=0; k<NoCommon; k++) CommonBonds[k] = 0;
                    
                    for (k=0; k<NoCommon; k++){
                        for (l=k+1; l<NoCommon; l++){
                            if (CNA_CheckBond(CommonNeighList[k]
                                              ,CommonNeighList[l])){
                                CommonBonds[k]++;
                                CommonBonds[l]++;
                                nobonds++;
                            }
                            
                        }
                    }
                    if (NoCommon==4 && nobonds==2){
                        for (k=0; k<NoCommon; k++){
                            maxbond = MAX(maxbond,CommonBonds[k]);
                            minbond = MIN(minbond,CommonBonds[k]);
                        }
                        if (minbond==1 && maxbond==1)nfcc++;
                        if (minbond==0 && maxbond==2)nhcp++;
                    }else if (NoCommon==5 && nobonds==5){
                        for (k=0; k<NoCommon; k++){
                            maxbond = MAX(maxbond,CommonBonds[k]);
                            minbond = MIN(minbond,CommonBonds[k]);
                        }
                        if (minbond==5 && maxbond==2)nico++;
                    }
                }
            }
            if (nfcc==12)pattern=CNA_FCC;
            else if (nfcc==6 && nhcp==6) pattern=CNA_HCP;
            else if (nico==12) pattern=CNA_ICOS;
            else pattern=CNA_OTHER;
        }
        else if (NoNeighs==14){
            int nbcc4,nbcc6;
            nbcc4=nbcc6=0;
            for (j=N_NonPairwise->idx[iatom];
                 j<N_NonPairwise->idx[iatom+1]; j++){
                NoCommon = 0;
                maxbond = 0;
                minbond = ATOM_COORDINATION_MAX;
                nobonds = 0;
                jatom = N_NonPairwise->list[j];
                for (k=N_NonPairwise->idx[iatom];
                     k<N_NonPairwise->idx[iatom+1]; k++){
                    katom = N_NonPairwise->list[k];
                    if(CNA_CheckBond(jatom,katom))
                        CommonNeighList[NoCommon++] = katom;
                }
                if (NoCommon==4 || NoCommon==6){
                    for (k=0; k<NoCommon; k++) CommonBonds[k] = 0;
                    
                    for (k=0; k<NoCommon; k++){
                        for (l=k+1; l<NoCommon; l++){
                            if (CNA_CheckBond(CommonNeighList[k]
                                              ,CommonNeighList[l])){
                                CommonBonds[k]++;
                                CommonBonds[l]++;
                                nobonds++;
                            }
                            
                        }
                    }
                    if (NoCommon == nobonds){
                        for (k=0; k<NoCommon; k++){
                            maxbond = MAX(maxbond,CommonBonds[k]);
                            minbond = MIN(minbond,CommonBonds[k]);
                        }
                        if (maxbond==2 && minbond==2){
                            if(NoCommon==4) nbcc4++;
                            if(NoCommon==6) nbcc6++;
                        }
                    }
                }
            }
            if (nbcc4==6 && nbcc6==8)pattern=CNA_BCC;
            else pattern=CNA_OTHER;
        }

        CNA_pattern[iatom] = pattern;
        CNA_hist[pattern]++;
    }
    ComputeCNA = 1;
    Neighborlist_free_nonpairwise_image
    (N_NonPairwise);
    return;
} /* end AppendCNA() */

/* Check wether two atoms are bonded or not*/
bool CNA_CheckBond(int iatom,int jatom){
    int i;
    for (i=N->idx[iatom];
         i<N->idx[iatom+1]; i++)
        if (N->list[i]==jatom) return (TRUE);
    for (i=N->idx[jatom];
         i<N->idx[jatom+1]; i++)
        if (N->list[i]==iatom) return (TRUE);
    return (FALSE);
}/* end CNA_CheckBond() */

bool assign_CNA_color (int iw)
{
    register int i;
    if (!n[iw].xtal_mode)
    {
        printf ("CNA color works only in Xtal mode.\n");
        return (FALSE);
    }
    strcpy(AX_title[iw], str4(fbasename, " (", "CNA", ")"));
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    for (i=0; i<np; i++)
    {
        AX_3D_AssignRGB(B->BALL[i], ATOM_CNA_COLOR[CNA_pattern[i]].r,
                        ATOM_CNA_COLOR[CNA_pattern[i]].g,
                        ATOM_CNA_COLOR[CNA_pattern[i]].b);
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]])
        * n[iw].atom_r_ratio;
    }
    n[iw].color_mode = COLOR_MODE_CNA;
    if (n[iw].bond_mode)
    {
        bond_xtal_origin_update(iw);
        bond_atom_color_update(iw);
    }
    else
    {
        n[iw].bond_xtal_origin_need_update = TRUE;
        n[iw].bond_atom_color_need_update = TRUE;
    }
    return (TRUE);
}

void print_CNA_histogram()
{
    int i;
    char *names[5];
    names[CNA_FCC]=" FCC   ";
    names[CNA_BCC]=" BCC   ";
    names[CNA_HCP]=" HCP   ";
    names[CNA_ICOS]=" ICOS  ";
    names[CNA_OTHER]=" Other ";
    printf ("----------- Common Neghbor Analysis Statistics ------------\n");
    printf ("Pattern.     Count  Percentage    R     G     B    Name\n");
    for (i=0; i<5; i++)
        if (CNA_hist[i] > 0)
            printf("%s %10d   %6.2f%%    %.3f %.3f %.3f  %s\n",
                   names[i], CNA_hist[i], 100.*CNA_hist[i]/np,
                   ATOM_CNA_COLOR[i].r, ATOM_CNA_COLOR[i].g,
                   ATOM_CNA_COLOR[i].b,
                   ATOM_CNA_COLOR[i].name );
    printf ("-----------------------------------------------------------\n");
    return;
}

