#include<iostream.h>
#include<math.h>
#include<stdlib.h>
#include<fstream.h>
#include<iomanip.h>

#include<windows.h>
#include<gl\gl.h>
#include<gl\glu.h>
#include<gl\glaux.h>
#include<gl\glut.h>

//some constants
double AVO=6.0225e23;
double BOLZ=1.38054;
double PI=3.14159265;
//read H0 matrix.
int VYON=0;  //with or without VEL
double H0[3][3];
int snhg=0;
int BigB=0;
int Line=0;

//show configuration or not.
int TryShow=0;

double Y;
int Nx,Ny,Nz,Vac,Def;
int Tot_num,Tot_num1;

double *X_mol;//global value.
double *Y_mol;
double *Z_mol;

double *XV_mol;
double *YV_mol;
double *ZV_mol;

double *XSE_mol;
double *YSE_mol;
double *ZSE_mol;

double *XTH_mol;
double *YTH_mol;	
double *ZTH_mol;

double *XFO_mol;
double *YFO_mol;
double *ZFO_mol;

double *XFI_mol;
double *YFI_mol;
double *ZFI_mol;

double *ForX_mol;
double *ForY_mol;	
double *ForZ_mol;

double *DISP_X;
double *DISP_Y;
double *DISP_Z;

double *XOLD_mol;
double *YOLD_mol;
double *ZOLD_mol;


int *IDIST;
int *NABORS;
int *LIST;

int ISUM;
double XSUM;
double SUME;
double SUMV;
double SUMVEL;
int NABTOT;
double TOTE;
double TOTV;
double EK;
double ESHFT;
double ESHFTA;
double FSHFT;
double RC;
double RDMAX;
double RMAX;
double RRMAX;
double RRMAX6;


int IFLG;
int KB;
int KSAVE;
int KSORT;
int KWRITE;
int MAXKB;
double TDIST;
double XDIST;
double TSTEP;


double TR;
double DR;
double DELTA;
double SIGMA;
double EPSI;
double WIMOL;
double RDEL;
double T;
double VO;
double VSCAL;
double VOLCC;
double RLIST;
int	NY;
double CORE;
double CORV;
double DE;
double DP;
double VOL;
double VSCALE;
int NMAX;


double RTR;//The defination is not clear.
double AHEAT;
double F0,F1,F2,F3,F4,F5;


double FKB;
double TMP;
double ENR;
double VIR;
double PRES;
double E1;
double ETOT;
double P1;
double RLTIM;

double ABSDIS=0.00001;
double ABSDIS1=0.0001;

/*                                                    the parameters in OpenGl                                                           */

const GLfloat lightPosition[] = {10.0,10.0,10.0,0.0};   
const GLfloat whiteLight[] = {0.8,0.8,0.8,1.0}; 
GLfloat matSpecular [] = {0.3,0.3,0.3,1.0};   
GLfloat matShininess [] = {20.0};   
GLfloat matEmission []={0,0,1.0,1.0} ;//= {238/256.0,199/256.0,16/256.0,1.0};Yellow   
GLfloat AMP_XYZ;
GLfloat	HUGEXTZ=20;
GLfloat YFourth=Y/4;

GLfloat	xrot=0;								
GLfloat	yrot=0;	
GLfloat	xtrf=0;								
GLfloat	ytrf=0;							
						
GLfloat xspeed=0.4f;								
GLfloat yspeed=0.4f;														
GLfloat xtranf=0.4f;
GLfloat ytranf=0.4f;
int Num_p=0;
int Num_pn=0;

FCC(double Y,int Nx,int Ny,int Nz,int vacancy,int defect,int Tot_num)//rectangle Nx*Ny*Nz
{
	double Z=Y/2;
	int P=0;
	int PP=0;
	int PPP=0;
	int XC,YC,ZC,Diff;
	double XCC,YCC,ZCC;

	*X_mol=0;*Y_mol=0;*Z_mol=0;
	*(X_mol+1)=Z;*(Y_mol+1)=0;*(Z_mol+1)=Z;
	*(X_mol+2)=Z;*(Y_mol+2)=Z;*(Z_mol+2)=0;
	*(X_mol+3)=0;*(Y_mol+3)=Z;*(Z_mol+3)=Z;
	
	for(int nz=0;nz<Nz;nz=nz+1)
	{
		for(int p=0;p<4;p=p+1)
		{
			*(X_mol+p+P)=*(X_mol+p);
			*(Y_mol+p+P)=*(Y_mol+p);
			*(Z_mol+p+P)=*(Z_mol+p)+Y*(nz);
		}
		P=P+4;
	}

	for(int ny=0;ny<Ny;ny=ny+1)
	{
		for(int pp=0;pp<4*Nz;pp=pp+1)
		{
			*(X_mol+pp+PP)=*(X_mol+pp);
			*(Y_mol+pp+PP)=*(Y_mol+pp)+Y*(ny);
			*(Z_mol+pp+PP)=*(Z_mol+pp);
		}
		PP=PP+4*Nz;
	}

	for(int nx=0;nx<Nx;nx=nx+1)
	{
		for(int ppp=0;ppp<4*Nz*Ny;ppp=ppp+1)
		{
			*(X_mol+ppp+PPP)=*(X_mol+ppp)+Y*(nx);
			*(Y_mol+ppp+PPP)=*(Y_mol+ppp);
			*(Z_mol+ppp+PPP)=*(Z_mol+ppp);
		}
		PPP=PPP+4*Nz*Ny;
	}

	Diff=vacancy-defect;//point defect
	if(Diff>0)
	{
		for(int lef=0;lef<Diff;lef=lef+1)
		{
			XC=(rand()/32767.0)*(4*Nx*Ny*Nz-1-lef);//The cost of equal probability.
			for(int i=XC;i<4*Nx*Ny*Nz-lef;i=i+1)
				*(X_mol+i)=*(X_mol+i+1);
		}
	}
	if(Diff<0)
	{
		for(int myw=1;myw<=abs(Diff);myw=myw+1)
		{
			 XC=(rand()/32767.0)*Nx;
			 YC=(rand()/32767.0)*Ny;
	    	 ZC=(rand()/32767.0)*Nz;
			 XCC=(rand()/32767.0)*Y;
			 YCC=(rand()/32767.0)*Y;
			 ZCC=(rand()/32767.0)*Y;
			 *(X_mol+myw+4*Nx*Ny*Nz)=XCC+Y*(XC);
			 *(Y_mol+myw+4*Nx*Ny*Nz)=YCC+Y*(YC);
			 *(Z_mol+myw+4*Nx*Ny*Nz)=ZCC+Y*(ZC);
		}
	}
}

int READCFGINF(char name[])
{
	char MYW[200];
	int sn=0;
	int Number=0;
	double Ans=1;
	int vyon=1;
	ifstream CFG (name);
	if (! CFG.is_open())
	{ 
		cout << "Error opening file"<<endl;
		return 1; 
	}
	while (! CFG.eof() ) 
	{
		sn=sn+1;
		CFG.getline (MYW,150);
		if(MYW[0]!='#')
		{
			if(sn==1)//Define this file and read the number.
			{
		    	if(MYW[0]=='N'&&MYW[1]=='u'&&MYW[2]=='m')
			    	cout<<"This is a simplified CFG file."<<endl;
				int kb=22;
				int pm=0;
	     		while(MYW[kb]!='\0')
				{
			    	pm=pm+1;
			    	kb=kb+1;
				}
		    	for(int l=0;l<pm;l=l+1)
			    	Number=(int(MYW[22+l])-48)*pow(10,pm-l-1)+Number;
				Tot_num=Number;
	    		cout<<"The total Number is "<<Tot_num<<endl;
			}


			if(MYW[0]=='A'&&MYW[1]==' '&&MYW[2]=='=')//read user's set.
			{
				Ans=0;int kb1=4;int pm1=0;int pm2=0;
            	while(MYW[kb1]!=' ')
				{
					pm1=pm1+1;
	            	kb1=kb1+1;
				}
            	for(int kk=4;kk<kb1;kk=kk+1)
				{
            		if(MYW[kk]=='.')
		        	pm2=kk;
				}
            	if(pm2==0)
				{
             		for(int kkk=0;kkk<pm1;kkk=kkk+1)
		            	Ans=Ans+(int(MYW[4+kkk])-48)*pow(10.0,pm1-kkk-1);
				}
            	if(pm2!=0)
				{
            		pm2=pm2-4;
	            	if(pm2==0)
					{
             			for(int kkk2=1;kkk2<pm1;kkk2=kkk2+1)
			        	Ans=Ans+(int(MYW[4+kkk2])-48)*pow(10.0,-kkk2);
					}
					if(pm2>0&&pm2<(pm1-1))
					{
						for(int kkk3=0;kkk3<pm2;kkk3=kkk3+1)
							Ans=Ans+(int(MYW[kkk3+4])-48)*pow(10.0,pm2-1-kkk3);
	             		for(int kkk4=pm2;kkk4<(pm1-1);kkk4=kkk4+1)
		            		Ans=Ans+(int(MYW[kkk4+5])-48)*pow(10.0,pm2-kkk4-1);
					}
					if(pm2==(pm1-1))
					{
						for(int kkk5=0;kkk5<pm2;kkk5=kkk5+1)
							Ans=Ans+(int(MYW[kkk5+4])-48)*pow(10.0,pm2-1-kkk5);
					}
				}				
				cout<<"The scale you set is "<<Ans<<" Angstrom."<<endl;
				SIGMA=Ans;
			}
			if(MYW[0]=='H'&&MYW[1]=='0'&&MYW[2]=='(')
			{
				int n1=(int(MYW[3])-48-1);
				int n2=(int(MYW[5])-48-1);

				double Ans1=0;int kb1=10;int pm1=0;int pm2=0;
            	while(MYW[kb1]!=' ')
				{
					pm1=pm1+1;
	            	kb1=kb1+1;
				}
            	for(int kk=10;kk<kb1;kk=kk+1)
				{
            		if(MYW[kk]=='.')
		        	pm2=kk;
				}
            	if(pm2==0)
				{
             		for(int kkk=0;kkk<pm1;kkk=kkk+1)
		            	Ans1=Ans1+(int(MYW[10+kkk])-48)*pow(10.0,pm1-kkk-1);
				}
            	if(pm2!=0)
				{
            		pm2=pm2-10;
	            	if(pm2==0)
					{
             			for(int kkk2=1;kkk2<pm1;kkk2=kkk2+1)
			        	Ans1=Ans1+(int(MYW[10+kkk2])-48)*pow(10.0,-kkk2);
					}
					if(pm2>0&&pm2<(pm1-1))
					{
						for(int kkk3=0;kkk3<pm2;kkk3=kkk3+1)
							Ans1=Ans1+(int(MYW[kkk3+10])-48)*pow(10.0,pm2-1-kkk3);
	             		for(int kkk4=pm2;kkk4<(pm1-1);kkk4=kkk4+1)
		            		Ans1=Ans1+(int(MYW[kkk4+10+1])-48)*pow(10.0,pm2-kkk4-1);
					}
					if(pm2==(pm1-1))
					{
						for(int kkk5=0;kkk5<pm2;kkk5=kkk5+1)
							Ans1=Ans1+(int(MYW[kkk5+10])-48)*pow(10.0,pm2-1-kkk5);
					}
				}	
				H0[n1][n2]=Ans1;
				cout<<"     H0("<<n1+1<<","<<n2+1<<") = "<<H0[n1][n2]<<" Angstrom."<<endl;
			}

			if(MYW[0]=='.'&&MYW[1]=='N'&&MYW[2]=='O')
			{
				vyon=0;
				if(vyon==0)
					cout<<"There is no velocity in the file."<<endl;
				else
					cout<<"There is velocity in the file."<<endl;
				VYON=vyon;
			}
			if(MYW[0]=='e'&&MYW[1]=='n'&&MYW[2]=='t')
			{
				int kb1=14;int pm1=0;double Ans2=0;
				while(MYW[kb1]!='\0')
				{
					pm1=pm1+1;
	            	kb1=kb1+1;
				}
				for(int kkk=0;kkk<pm1;kkk=kkk+1)
					Ans2=Ans2+(int(MYW[14+kkk])-48)*pow(10.0,pm1-kkk-1);
				cout<<"The entry_count is "<<Ans2<<endl;
			}
		}
     }
	 CFG.close();
     return 0;
}

int READFIL(char name[],int VYON)
{
	char MYW[200];
	ifstream CFG (name);
	if (! CFG.is_open())
	{ 
		cout << "Error opening file"<<endl;
		return 1; 
	}
	while (! CFG.eof() ) 
	{
		CFG.getline (MYW,180);


		if(snhg==1000)
			BigB=1000;
		if(MYW[2]=='\0')//find the critical line.
			snhg=1000;


		if(BigB==1000&&MYW[0]!='\0')        //read the data from file.
		{
		//	cout<<BigB<<endl;
			int ABC=3;
			if(VYON==1)
				ABC=6;
			int kb1=0;
			int pm1=0;
			int pm2=-1;
			double AnsFil=0;
			int ILikeM=0;
			int kbtmp=0;
			while(MYW[kb1]!='\0')
			{
			    kbtmp=kb1;		
				pm1=0;
				AnsFil=0;
				while((MYW[kb1]!=' ')&&(MYW[kb1]!='\0'))
				{
					pm1=pm1+1;
	            	kb1=kb1+1;
				}
//				cout<<kb1<<","<<pm1<<",";
				if(MYW[kb1]==' '||MYW[kb1]=='\0')
					ILikeM=ILikeM+1;
//				cout<<ILikeM<<","<<Line<<",";
            	for(int kk=kbtmp;kk<kb1;kk=kk+1)
				{
            		if(MYW[kk]=='.')
		        	pm2=kk;
				}
            	if(pm2==-1)
				{
             		for(int kkk=0;kkk<pm1;kkk=kkk+1)
		            	AnsFil=AnsFil+(int(MYW[kbtmp+kkk])-48)*pow(10.0,pm1-kkk-1);
				}
            	if(pm2!=-1)
				{
            		pm2=pm2-kbtmp;
	            	if(pm2==0)
					{
             			for(int kkk2=1;kkk2<pm1;kkk2=kkk2+1)
			        	AnsFil=AnsFil+(int(MYW[kbtmp+kkk2])-48)*pow(10.0,-kkk2);
					}
					if(pm2>0&&pm2<(pm1-1))
					{
						for(int kkk3=0;kkk3<pm2;kkk3=kkk3+1)
							AnsFil=AnsFil+(int(MYW[kkk3+kbtmp])-48)*pow(10.0,pm2-1-kkk3);
	             		for(int kkk4=pm2;kkk4<(pm1-1);kkk4=kkk4+1)
		            		AnsFil=AnsFil+(int(MYW[kkk4+kbtmp+1])-48)*pow(10.0,pm2-kkk4-1);
					}
					if(pm2==(pm1-1))
					{
						for(int kkk5=0;kkk5<pm2;kkk5=kkk5+1)
							AnsFil=AnsFil+(int(MYW[kkk5+kbtmp])-48)*pow(10.0,pm2-1-kkk5);
					}
				}
//				cout<<AnsFil<<endl;
				if(ABC==6)
				{
					if(ILikeM%ABC==1)
				    	*(X_mol+Line)=AnsFil;
			    	if(ILikeM%ABC==2)
				    	*(Y_mol+Line)=AnsFil;
				    if(ILikeM%ABC==3)
				    	*(Z_mol+Line)=AnsFil;
			    	if(ILikeM%ABC==3)
	     		    	*(XV_mol+Line)=AnsFil;
			     	if(ILikeM%ABC==5)
	     		    	*(YV_mol+Line)=AnsFil;
			    	if(ILikeM%ABC==0)
	     		    	*(ZV_mol+Line)=AnsFil;
				}
				if(ABC==3)
				{
					if(ILikeM%ABC==1)
				    	*(X_mol+Line)=AnsFil;
			    	if(ILikeM%ABC==2)
				    	*(Y_mol+Line)=AnsFil;
				    if(ILikeM%ABC==0)
				    	*(Z_mol+Line)=AnsFil;
				}
				if(MYW[kb1]!='\0')
					kb1=kb1+1;
			}

			Line=Line+1;
		}
	}
	CFG.close();
	return 0;
}


int WRITECFG(char name[],int Tot_num,double SIGMA,int Nx,int Ny,int Nz,double Y,int FLGVEL)
{
	ofstream out(name);
	if(!out)
	{
		cout<<"Cannot open!"<<endl;
		return 1;
	}
	out<<"Number of particles = "<<Tot_num<<endl;
	out<<"A = "<<SIGMA<<" Angstrom"<<endl;
	out<<"H0(1,1) = "<<Nx*Y<<" A"<<endl;
	out<<"H0(1,2) = "<<0<<" A"<<endl;
	out<<"H0(1,3) = "<<0<<" A"<<endl;

	out<<"H0(2,1) = "<<0<<" A"<<endl;
	out<<"H0(2,2) = "<<Ny*Y<<" A"<<endl;
	out<<"H0(2,3) = "<<0<<" A"<<endl;

	out<<"H0(3,1) = "<<0<<" A"<<endl;
	out<<"H0(3,2) = "<<0<<" A"<<endl;
	out<<"H0(3,3) = "<<Nz*Y<<" A"<<endl;

	if(FLGVEL==0)
	{
		out<<".NO_VELOCITY."<<endl;
		out<<"entry_count = 3"<<endl;
		out<<"200.59"<<endl;
    	out<<"Hg"<<endl;

    	for(int kk=0;kk<Tot_num;kk=kk+1)
    		out<<*(X_mol+kk)<<" "<<*(Y_mol+kk)<<" "<<*(Z_mol+kk)<<endl;
	}
	if(FLGVEL==1)
	{
		out<<"entry_count = 6"<<endl;
		out<<"200.59"<<endl;
    	out<<"Hg"<<endl;

    	for(int kk=0;kk<Tot_num;kk=kk+1)
    		out<<*(X_mol+kk)<<" "<<*(Y_mol+kk)<<" "<<*(Z_mol+kk)<<*(XV_mol+kk)<<" "<<*(YV_mol+kk)<<" "<<*(ZV_mol+kk)<<endl;
	}
	out.close();
	return 0;
}

void BACKGROUND(void)
{
	glClearColor(1.0,1.0,1.0,1.0);
}

void MYDISPLAY()
{
   Num_p=0;
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);   
   glMatrixMode(GL_PROJECTION);   
   glLoadIdentity(); 
   glOrtho(-10,10,-10,10.0,-10,10.0);    
   glRotatef(xrot,1.0f,0.0f,0.0f);		
   glRotatef(yrot,0.0f,1.0f,0.0f);		
   glTranslatef(xtranf,0.0f,0.0f);
   glTranslatef(ytranf,0.0f,0.0f);

   YFourth=Y*0.6;
   glMatrixMode(GL_MODELVIEW);   
   glPushMatrix();   

   glEnable(GL_LINES);
   glBegin(GL_LINES);

   GLfloat npx=H0[0][0]*AMP_XYZ;
   glColor3f(1.0f,0,0);
   glVertex3f(0,0,0);
   glVertex3f(npx,0,0);

   GLfloat npy=H0[1][1]*AMP_XYZ;
   glColor3f(0,1.0f,0);
   glVertex3f(0,0,0);
   glVertex3f(0,npy,0);

   GLfloat npz=H0[2][2]*AMP_XYZ;
   glColor3f(0,0,1.0f);
   glVertex3f(0,0,0);
   glVertex3f(0,0,npz);

   glEnd();
   glFlush();

   glMaterialfv(GL_FRONT,GL_SPECULAR,matSpecular);   
   glMaterialfv(GL_FRONT,GL_SHININESS,matShininess);   
   glMaterialfv(GL_FRONT,GL_EMISSION,matEmission);   
   
   for(int myw=0;myw<Tot_num;myw=myw+1)
   {

	   GLfloat x1=(*(X_mol+myw))*AMP_XYZ/5;
	   GLfloat y1=(*(Y_mol+myw))*AMP_XYZ/5;
	   GLfloat z1=(*(Z_mol+myw))*AMP_XYZ/5;
	   glTranslatef(x1,y1,z1);
	   glColor3f(1.0,-1,-1);
	   glutSolidSphere(YFourth,16,16);
	   Num_p=Num_p+1;
	   glTranslatef(-x1,-y1,-z1);
   }
   if(Num_pn==0)
   {
	   cout<<"There are "<<Num_p<<" atoms have been painted."<<endl;
	   Num_pn=Num_pn+1;
   }
   glPopMatrix();   
   glFlush();
   xrot+=xspeed;						
   yrot+=yspeed;
   xtranf+=xtrf;
   ytranf+=ytrf;
   return;
}

void INIT()   
{   
    glClearColor(1.0,1.0,1.0,0);   
    glClearDepth(1.0);   
    glShadeModel(GL_SMOOTH);   
    glEnable(GL_LIGHTING);   
    glEnable(GL_LIGHT0);   
    glEnable(GL_DEPTH_TEST);   
    glMatrixMode(GL_MODELVIEW);   
    glLoadIdentity();   
    glLightfv(GL_LIGHT0,GL_POSITION,lightPosition);   
    glLightfv(GL_LIGHT0,GL_DIFFUSE,whiteLight);   
    glLightfv(GL_LIGHT0,GL_SPECULAR,whiteLight);   
}   


void MYRESHAPE(int w,int h)   
{   
    glViewport(0.0,0.0,(GLsizei) w,(GLsizei) h);
	HUGEXTZ=((w>h)?h:w)/800*20;
	double MYW=NMAX*Y;
	AMP_XYZ=HUGEXTZ/MYW;
}   

void KEYBOARD(unsigned char key,int x,int y)
{
	switch(key)
	{
	case 'a':
		xspeed=xspeed;	
		break;
	case 'd':
		xspeed=-xspeed;
		break;
	case 'w':
		yspeed=yspeed;
		break;
	case 's':
		yspeed=-yspeed;
		break;
	case 'p':
		AMP_XYZ=AMP_XYZ/1.11;
		break;
	case 'l':
		AMP_XYZ=1.11*AMP_XYZ;
		break;
	case 'h'://you can see atoms in xy plane or yz or xz plane by using 'h'\'j'
		xtranf=xtranf;
		break;
	case 'j':
		ytranf=ytranf;
		break;
	case 'b':
		YFourth=YFourth*2;
		break;
	case 'n':
		YFourth=YFourth/2;
		break;
	}
		glutPostRedisplay();
}

INIVEL(double RTR,double AHEAT,int Tot_num)
{
	double SUMX=0;
	double SUMY=0;
	double SUMZ=0;
	double KIN=0;//kinetic energe
	for(int q=0;q<Tot_num;q=q+1)
	{
		double XV=(rand()/32767.0-0.5)*2;
		double YV=(rand()/32767.0-0.5)*2;
		double ZV=(rand()/32767.0-0.5)*2;
		double Total=sqrt(XV*XV+YV*YV+ZV*ZV);
		*(XV_mol+q)=XV*RTR/Total;
		*(YV_mol+q)=YV*RTR/Total;
		*(ZV_mol+q)=ZV*RTR/Total;
		SUMX=*(XV_mol+q)+SUMX;
		SUMY=*(YV_mol+q)+SUMY;
		SUMZ=*(ZV_mol+q)+SUMZ;
	}
	for(int r=0;r<Tot_num;r=r+1)//momenta conservation
	{
		*(XV_mol+r)=*(XV_mol+r)-SUMX/Tot_num;
		*(YV_mol+r)=*(YV_mol+r)-SUMY/Tot_num;
		*(ZV_mol+r)=*(ZV_mol+r)-SUMZ/Tot_num;
		KIN=KIN+(*(XV_mol+r))*(*(XV_mol+r))+(*(YV_mol+r))*(*(YV_mol+r))+(*(ZV_mol+r))*(*(ZV_mol+r));
	}
	for(int s=0;s<Tot_num;s=s+1)//desired temperature
	{
		*(XV_mol+s)=*(XV_mol+s)*sqrt(AHEAT/KIN);
		*(YV_mol+s)=*(YV_mol+s)*sqrt(AHEAT/KIN);
		*(ZV_mol+s)=*(ZV_mol+s)*sqrt(AHEAT/KIN);
	}
}


INIACC(double DELTA,int Tot_num)
{
	for(int myw=0;myw<Tot_num;myw=myw+1)
	{
		*(XSE_mol+myw)=*(ForX_mol+myw)*0.5*DELTA*DELTA;
		*(YSE_mol+myw)=*(ForY_mol+myw)*0.5*DELTA*DELTA;
		*(ZSE_mol+myw)=*(ForZ_mol+myw)*0.5*DELTA*DELTA;
		*(XOLD_mol+myw)=*(X_mol+myw);
		*(YOLD_mol+myw)=*(Y_mol+myw);
		*(ZOLD_mol+myw)=*(Z_mol+myw);
	}
}


PREDCT(int Tot_num)//predict the new configuration
{
	for(int t=0;t<Tot_num;t=t+1)
	{
		*(X_mol+t)=*(X_mol+t)+*(XV_mol+t)+*(XSE_mol+t)+*(XTH_mol+t)+*(XFO_mol+t)+*(XFI_mol+t);
		*(Y_mol+t)=*(Y_mol+t)+*(YV_mol+t)+*(YSE_mol+t)+*(YTH_mol+t)+*(YFO_mol+t)+*(YFI_mol+t);
		*(Z_mol+t)=*(Z_mol+t)+*(ZV_mol+t)+*(ZSE_mol+t)+*(ZTH_mol+t)+*(ZFO_mol+t)+*(ZFI_mol+t);

		*(XV_mol+t)=*(XV_mol+t)+2*(*(XSE_mol+t))+3*(*(XTH_mol+t))+4*(*(XFO_mol+t))+5*(*(XFI_mol+t));
		*(YV_mol+t)=*(YV_mol+t)+2*(*(YSE_mol+t))+3*(*(YTH_mol+t))+4*(*(YFO_mol+t))+5*(*(YFI_mol+t));
		*(ZV_mol+t)=*(ZV_mol+t)+2*(*(ZSE_mol+t))+3*(*(ZTH_mol+t))+4*(*(ZFO_mol+t))+5*(*(ZFI_mol+t));

		*(XSE_mol+t)=*(XSE_mol+t)+3*(*(XTH_mol+t))+6*(*(XFO_mol+t))+10*(*(XFI_mol+t));
		*(YSE_mol+t)=*(YSE_mol+t)+3*(*(YTH_mol+t))+6*(*(YFO_mol+t))+10*(*(YFI_mol+t));
		*(ZSE_mol+t)=*(ZSE_mol+t)+3*(*(ZTH_mol+t))+6*(*(ZFO_mol+t))+10*(*(ZFI_mol+t));

		*(XTH_mol+t)=*(XTH_mol+t)+4*(*(XFO_mol+t))+10*(*(XFI_mol+t));
		*(YTH_mol+t)=*(YTH_mol+t)+4*(*(YFO_mol+t))+10*(*(YFI_mol+t));
		*(ZTH_mol+t)=*(ZTH_mol+t)+4*(*(ZFO_mol+t))+10*(*(ZFI_mol+t));

		*(XFO_mol+t)=*(XFO_mol+t)+5*(*(XFI_mol+t));
		*(YFO_mol+t)=*(YFO_mol+t)+5*(*(YFI_mol+t));
		*(ZFO_mol+t)=*(ZFO_mol+t)+5*(*(ZFI_mol+t));

		*(ForX_mol+t)=0;
		*(ForY_mol+t)=0;
		*(ForZ_mol+t)=0;
	}
}

CORR(double Y,int Nx,int Ny,int Nz,int Tot_num,double DELTA)//the correction step
{
	for(int w=0;w<Tot_num;w=w+1)
	{
		double XERR=*(XSE_mol+w)-0.5*DELTA*DELTA*(*(ForX_mol+w));
		double YERR=*(YSE_mol+w)-0.5*DELTA*DELTA*(*(ForY_mol+w));
		double ZERR=*(ZSE_mol+w)-0.5*DELTA*DELTA*(*(ForZ_mol+w));
		*(X_mol+w)=*(X_mol+w)-XERR*F0;
		*(XV_mol+w)=*(XV_mol+w)-XERR*F1;
		*(XSE_mol+w)=*(XSE_mol+w)-XERR*F2;
		*(XTH_mol+w)=*(XTH_mol+w)-XERR*F3;
		*(XFO_mol+w)=*(XFO_mol+w)-XERR*F4;
		*(XFI_mol+w)=*(XFI_mol+w)-XERR*F5;

		*(Y_mol+w)=*(Y_mol+w)-YERR*F0;
		*(YV_mol+w)=*(YV_mol+w)-YERR*F1;
		*(YSE_mol+w)=*(YSE_mol+w)-YERR*F2;
		*(YTH_mol+w)=*(YTH_mol+w)-YERR*F3;
		*(YFO_mol+w)=*(YFO_mol+w)-YERR*F4;
		*(YFI_mol+w)=*(YFI_mol+w)-YERR*F5;

		*(Z_mol+w)=*(Z_mol+w)-ZERR*F0;
		*(ZV_mol+w)=*(ZV_mol+w)-ZERR*F1;
		*(ZSE_mol+w)=*(ZSE_mol+w)-ZERR*F2;
		*(ZTH_mol+w)=*(ZTH_mol+w)-ZERR*F3;
		*(ZFO_mol+w)=*(ZFO_mol+w)-ZERR*F4;
		*(ZFI_mol+w)=*(ZFI_mol+w)-ZERR*F5;

		*(DISP_X+w)=*(DISP_X+w)-(*(X_mol+w))+(*(XOLD_mol+w));
		*(DISP_Y+w)=*(DISP_Y+w)-(*(Y_mol+w))+(*(YOLD_mol+w));
		*(DISP_Z+w)=*(DISP_Z+w)-(*(Z_mol+w))+(*(ZOLD_mol+w));

		double XX=*(X_mol+w);
		double YY=*(Y_mol+w);
		double ZZ=*(Z_mol+w);

		int XXnum=XX/(Y*Nx);//the check for the particle number conservation
		int YYnum=YY/(Y*Ny);
		int ZZnum=ZZ/(Y*Nz);
		XX=XX-XXnum*(Y*Nx);
		YY=YY-YYnum*(Y*Ny);
		ZZ=ZZ-ZZnum*(Y*Nz);
		if(XX<0)
			XX=XX+(Y*Nx);
		if(YY<0)
			YY=YY+(Y*Ny);
		if(ZZ<0)
			ZZ=ZZ+(Y*Nz);

		*(X_mol+w)=XX;
		*(Y_mol+w)=YY;
		*(Z_mol+w)=ZZ;
		*(XOLD_mol+w)=XX;
		*(YOLD_mol+w)=YY;
		*(ZOLD_mol+w)=ZZ;

	}
}

EVAL(double Y,int Nx,int Ny,int Nz,int Tot_num,int KB,int KSORT,double RDMAX,double RDEL,double RLIST,double RMAX,double TOTE,double TOTV,int NABTOT)
{
	int k=-1;
	int JBEGIN=0;
	int JEND=0;
	int LCHK=2;
	TOTE=0;
	TOTV=0;
	if(KB%KSORT==0)
		LCHK=1;
	for(int u=0;u<Tot_num-1;u=u+1)
	{
		if(LCHK==1)
		{
			*(NABORS+u)=k+1;
			JBEGIN=u+1;
			JEND=Tot_num-1;
		}
		if(LCHK==2)
		{
			JBEGIN=*(NABORS+u);
			JEND=*(NABORS+u+1)-1;
		}
		if(JBEGIN>JEND)
			return;
		double X1=*(X_mol+u);
		double Y1=*(Y_mol+u);
		double Z1=*(Z_mol+u);
		for(int v=JBEGIN;v<=JEND;v=v+1)
		{
			int v1=v;
			if(LCHK==2)
				v1=*(LIST+v);

			double XRR=X1-(*(X_mol+v1));
			double YRR=Y1-(*(Y_mol+v1));
			double ZRR=Z1-(*(Z_mol+v1));
			double XR=fabs(XRR);
			double YR=fabs(YRR);
			double ZR=fabs(ZRR);

			if(XR>(Y*Nx/2.0))//find the nearest image atom.
			{
				if(XRR>0)
					XR=XR-Y*Nx;
				else
					XR=XR+Y*Nx;
			}
			if(YR>(Y*Ny/2))
			{
				if(YRR>0)
					YR=YR-Y*Ny;
				else
					YR=YR+Y*Ny;
			}
			if(ZR>(Y*Nz/2))
			{
				if(ZRR>0)
					ZR=ZR-Y*Nz;
				else
					ZR=ZR+Y*Nz;
			}
			double RSQ=XR*XR+YR*YR+ZR*ZR;

			if(LCHK==1)
			{
				if(RSQ<=RDMAX)
				{
					int IJ=int(sqrt(RSQ)/RDEL+0.5);
					*(IDIST+IJ)=*(IDIST+IJ)+1;
					if(RSQ<=RLIST)
					{
						k=k+1;
						*(LIST+k)=v1;
					}
				}
			}
			if(RSQ<=RMAX)
			{
				double RSI=1.0/RSQ;
	     		double R6=pow(RSI,3);
	    		double RPL=48*R6*(R6-0.5)-sqrt(RSQ)*FSHFT;//Force momenta.
	    		double RP=RPL*RSI;

		    	double REPX=RP*XRR;
		    	*(ForX_mol+u)=*(ForX_mol+u)+REPX;
	    		*(ForX_mol+v1)=*(ForX_mol+v1)-REPX;
		    	double REPY=RP*YRR;
		    	*(ForY_mol+u)=*(ForY_mol+u)+REPY;
	    		*(ForY_mol+v1)=*(ForY_mol+v1)-REPY;
		    	double REPZ=RP*ZRR;
		    	*(ForZ_mol+u)=*(ForZ_mol+u)+REPZ;
	    		*(ForZ_mol+v1)=*(ForZ_mol+v1)-REPZ;

				TOTE=TOTE+4*R6*(R6-1)+ESHFT+sqrt(RSQ)*ESHFTA;//Why so?
		    	TOTV=TOTV-RPL;
			}
		}
		if(LCHK==1)
		{
			*(NABORS+Tot_num-1)=k+1;
			ISUM=ISUM+Tot_num/2.0;
			NABTOT=k+1;
		}
	}
}


INITIAL()//initial control parameters
{

	NABTOT=0;
    TOTE=0.0;
	TOTV=0.0;

    IFLG=0;
    KB=0;
    KSAVE=10;
    KSORT=10;
    KWRITE=10;
    MAXKB=2000;
    TDIST=0.0;
    XDIST=0.4;

    DELTA=0.001;
    SIGMA=3.405;
    EPSI=120;
    WIMOL=39.994;
	RC=2.5;
	RDEL=0.025;

	//calculation
	AHEAT=DELTA*DELTA*Tot_num*3*TR;
	RTR=DELTA*sqrt(3*TR);
	T=TR*EPSI;
	VOL=Tot_num/DR;
	VSCALE=1.0/(3*TR);
	NMAX=(H0[0][0]>H0[1][1])?H0[0][0]:H0[1][1];
	NMAX=(NMAX>H0[2][2])?NMAX:H0[2][2];
	RDMAX=(NMAX/2.0)*(NMAX/2.0);
	VOLCC=AVO*pow((SIGMA*1.0e-8),3)/DR;
	RLIST=(RC+0.3)*(RC+0.3);
	if(RLIST>RDMAX)
		RLIST=RDMAX;
	if(RC>(NMAX/2))
		RC=NMAX/2-0.35;
	RMAX=RC*RC;
	AMP_XYZ=20.0/(NMAX);

	//shifted-force constants
	RRMAX=1.0/sqrt(RMAX);
	RRMAX6=pow(RRMAX,6);
	ESHFT=RRMAX6*(28-52*RRMAX6);
	ESHFTA=48*RRMAX*RRMAX6*(RRMAX6-0.5);
	FSHFT=ESHFTA;
	
	F0=3/16.0;
	F1=251/360.0;
	F2=1.0;
	F3=11/18.0;
	F4=1/6.0;
	F5=1/60.0;

	TSTEP=sqrt(WIMOL*SIGMA*SIGMA/AVO/EPSI/BOLZ)*1.0e12;
	NY=int(NMAX/(2*RDEL)-1);

	//initialize sum accummulators
	ISUM=0;
    XSUM=0;
    SUME=0;
    SUMV=0;
	SUMVEL=0;
	TDIST=0;
	EK=0;
	//There is an array called IDIST[NY],but how can we get this array out this function.
   
	//correction for long-range interactions
	CORE=8*PI*DR*(1.0/(9.0*pow(RC,9))-1.0/(3.0*pow(RC,3)));
    CORV=98*PI*DR*(1.0/(6.0*pow(RC,3))-1.0/(9.0*pow(RC,9)));
    DE=CORE;
    DP=-VSCALE*CORV;
}


RDF(double DR,double RDEL,int NY)
{
	int YTot;
	double PDEN=DR*pow(RDEL,3);
	for(int zh=0;zh<NY;zh=zh+1)
	{
		if(*(IDIST+zh)!=0)
		{
			double V=(4*PI/3.0)*(3*zh*zh);
			int X=*(IDIST+zh);
			double GR=X/V/YTot/PDEN;
			double RRR=RDEL*(zh);
		}
		if(*(IDIST+zh)==0)
			return;
	}
}

EQBRAT(int KB,int IFLG,double AHEAT,double SUMVEL,int TDIST,int XDIST,int Tot_num,int NY)
{
	double HEAT;
	if(TDIST<=XDIST)
	{
		HEAT=sqrt(AHEAT/SUMVEL);
		for(int yw=0;yw<Tot_num;yw=yw+1)
		{
			*(XV_mol+yw)=(*(XV_mol+yw))*HEAT;
			*(YV_mol+yw)=(*(YV_mol+yw))*HEAT;
			*(ZV_mol+yw)=(*(ZV_mol+yw))*HEAT;
		}
	}
	if(TDIST>XDIST)
	{
		if(IFLG>=0)
		{
			IFLG=-1;
			KB=0;
			ISUM=0;
			SUME=0;
			SUMV=0;
			XSUM=0;
			for(int mw=0;mw<NY;mw=mw+1)
				*(IDIST+mw)=0;
		}
		if(KB>=100)
		{
			IFLG=1;
		}
	}
}


int main(int argc,char *argv[])
{
	TR=1.1;
    DR=4;

	/*                                              Initialize configuration velocity and accerlation                                         */
    
	char name[100];
	int IPOOP=0;  //input or output
	cout<<"Would you want to input the configuration from file,1 respect for Yes,0 is No?"<<endl;
	cin>>IPOOP;
	if(IPOOP==0)
	{
		cout<<"input Nx,Ny,Nz,vacancy,defect"<<endl;
		cin>>Nx>>Ny>>Nz>>Vac>>Def;
		Tot_num=4*Nx*Ny*Nz-Vac+Def;
		Tot_num1=4*Nx*Ny*Nz+Def+1;//new memeory to ensure the size is bigger than Tot_num.It's the smallest one.*/
	    Y=pow(Tot_num/(DR*Nx*Ny*Nz),1.0/3);
		H0[0][0]=Nx*Y;H0[0][1]=0;H0[0][2]=0;
	    H0[1][0]=0;H0[1][1]=Ny*Y;H0[1][2]=0;
    	H0[2][0]=0;H0[2][1]=0;H0[2][2]=Nz*Y;
	}
	if(IPOOP==1)
	{
		cout<<"Input the name of file."<<endl;
		cin>>name;
		READCFGINF(name);//Input the information of the file£ºgeuss Tot_num;geuss VYON;
		Tot_num1=Tot_num+2;
	}


/*                                                     Got memory according to Tot_num1                                             */

	
	    X_mol=new double[Tot_num1];
    	Y_mol=new double[Tot_num1];
    	Z_mol=new double[Tot_num1];

    	XV_mol=new double[Tot_num1];
     	YV_mol=new double[Tot_num1];
    	ZV_mol=new double[Tot_num1];

        XSE_mol=new double[Tot_num1];
    	YSE_mol=new double[Tot_num1];
    	ZSE_mol=new double[Tot_num1];
	
    	XTH_mol=new double[Tot_num1];
    	YTH_mol=new double[Tot_num1];	
    	ZTH_mol=new double[Tot_num1];
	
    	XFO_mol=new double[Tot_num1];
    	YFO_mol=new double[Tot_num1];
    	ZFO_mol=new double[Tot_num1];

    	XFI_mol=new double[Tot_num1];
    	YFI_mol=new double[Tot_num1];
    	ZFI_mol=new double[Tot_num1];	

    	ForX_mol=new double[Tot_num1];
    	ForY_mol=new double[Tot_num1];	
    	ForZ_mol=new double[Tot_num1];

        DISP_X=new double[Tot_num1];
        DISP_Y=new double[Tot_num1];
        DISP_Z=new double[Tot_num1];

        XOLD_mol=new double[Tot_num1];
        YOLD_mol=new double[Tot_num1];
        ZOLD_mol=new double[Tot_num1];

    	NABORS=new int[Tot_num1];
    	IDIST=new int[NY];
    	int p=Tot_num*50;
        LIST=new int[p];

    	for(int pp=0;pp<Tot_num1;pp++)
		{
    		*(X_mol+pp)=0;
    		*(Y_mol+pp)=0;
    		*(Z_mol+pp)=0;
  
    		*(XV_mol+pp)=0;
    		*(YV_mol+pp)=0;
    		*(ZV_mol+pp)=0;

    		*(XSE_mol+pp)=0;
     		*(YSE_mol+pp)=0;
    		*(ZSE_mol+pp)=0;

    		*(XTH_mol+pp)=0;
    		*(YTH_mol+pp)=0;
    		*(ZTH_mol+pp)=0;

     		*(XFO_mol+pp)=0;
    		*(YFO_mol+pp)=0;
    		*(ZFO_mol+pp)=0;

    		*(XFI_mol+pp)=0;
    		*(YFI_mol+pp)=0;
	    	*(ZFI_mol+pp)=0;

        	*(ForX_mol+pp)=0;
        	*(ForY_mol+pp)=0;
        	*(ForZ_mol+pp)=0;
 
            *(DISP_X+pp)=0;
            *(DISP_Y+pp)=0;
            *(DISP_Z+pp)=0;
	
    		*(XOLD_mol+pp)=0;
    		*(YOLD_mol+pp)=0;
    		*(ZOLD_mol+pp)=0;


        	*(NABORS+pp)=0;
		}
    	for(int pp1=0;pp1<NY;pp1++)
    		*(IDIST+pp1)=0;
    	for(int pp2=0;pp2<50*Tot_num;pp2++)
    		*(LIST+pp2)=0;


/*                                                    Enter into the option                                              */

		if(IPOOP==0)
		{
			INITIAL();//Initialize all numbers
			FCC(Y,Nx,Ny,Nz,Vac,Def,Tot_num);
		}
		if(IPOOP==1)
		{
			READFIL(name,VYON);

			double SIGMAtmp=SIGMA;
			int num_Top=0;
			int NYZ=0;
			int NXZ=0;
			int NXY=0;
			int Num_Fir=0;
			for(int zj=0;zj<Tot_num;zj=zj+1)
			{
				while((fabs(*(X_mol+zj))>ABSDIS||fabs(*(Y_mol+zj))>ABSDIS||fabs(*(Z_mol+zj))>ABSDIS)&&Num_Fir==0)
				{
					num_Top=zj;
					Num_Fir=Num_Fir+1;
				}
			}
			for(int zj1=0;zj1<Tot_num;zj1=zj1+1)
			{
				if(fabs(*(X_mol+num_Top)-(*(X_mol+zj1)))<ABSDIS1)
					NYZ=NYZ+1;
				if(fabs(*(Y_mol+num_Top)-(*(Y_mol+zj1)))<ABSDIS1)
					NXZ=NXZ+1;
				if(fabs(*(Z_mol+num_Top)-(*(Z_mol+zj1)))<ABSDIS1)
					NXY=NXY+1;
			}
//			cout<<NYZ<<","<<NXZ<<","<<NXY<<endl;			
			Nx=int(double(Tot_num)/NYZ+0.2);
			Ny=int(double(Tot_num)/NXZ+0.2);
			Nz=int(double(Tot_num)/NXY+0.2);
			double YXYZ[3];
			YXYZ[0]=H0[0][0]/Nx;
			YXYZ[1]=H0[1][1]/Ny;
			YXYZ[2]=H0[2][2]/Nz;
			double H0MAX=YXYZ[0];
			double H0MIN=YXYZ[0];
			double H0MID=YXYZ[0];
			for(int iii=1;iii<=2;iii=iii+1)
			{
				H0MAX=(H0MAX>YXYZ[iii])?H0MAX:YXYZ[iii];
				H0MIN=(H0MIN<YXYZ[iii])?H0MIN:YXYZ[iii];
			}
			H0MID=YXYZ[0]+YXYZ[1]+YXYZ[2]-H0MAX-H0MIN;
			Y=H0MID;

//			cout<<Nx<<","<<Ny<<","<<Nz<<","<<Y<<endl;
			INITIAL();
			SIGMA=SIGMAtmp;
		}
		if(IPOOP==0||(IPOOP==1&VYON==0))
			INIVEL(RTR,AHEAT,Tot_num);
		if(IPOOP==0)
		{
			cout<<"The data is being writen to hardisk......"<<endl;
			cout<<"Rename the new file."<<endl;
			cin>>name;
			cout<<"Would you like to output the velocity into cfg?1 is respect for Yes,0 is No."<<endl;
			cin>>VYON;
			WRITECFG(name,Tot_num,SIGMA,Nx,Ny,Nz,Y,VYON);
		}

	cout<<"Would you like to show the atoms? 1 is respect for Yes.0 is No."<<endl;
	cin>>TryShow;
	if(TryShow==1)
	{
		glutInit(&argc,argv);
    	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGBA|GLUT_DEPTH);   
    	glutInitWindowSize(800,800);   
	   	glutInitWindowPosition(150,150);   
   		glutCreateWindow("FCC");
	   	BACKGROUND();
        glutDisplayFunc(MYDISPLAY);
        glutReshapeFunc(MYRESHAPE);  
        glutKeyboardFunc(KEYBOARD);
       	INIT();
       	glutMainLoop();
	}
	EVAL(Y,Nx,Ny,Nz,Tot_num,KB,KSORT,RDMAX,RDEL,RLIST,RMAX,TOTE,TOTV,NABTOT);//for vacancy and defect or form,this function is not correct.Please think this tommorrow. 
    INIACC(DELTA,Tot_num);

/*	double RMYWA=0;
	double RTY=0;
	for(int mm=1;mm<400;mm++)
	{
		KB=KB+1;
		cout<<KB<<":"<<endl;
		PREDCT(Tot_num);
		EVAL(Y,Nx,Ny,Nz,Tot_num,KB,KSORT,RDMAX,RDEL,RLIST,RMAX,TOTE,TOTV,NABTOT);
		CORR(Y,Nx,Ny,Nz,Tot_num,DELTA);
		double RMYW=(*(X_mol+0)-(*(X_mol+1)))*(*(X_mol+0)-(*(X_mol+1)))+(*(Y_mol+0)-(*(Y_mol+1)))*(*(Y_mol+0)-(*(Y_mol+1)))+(*(Z_mol+0)-(*(Z_mol+1)))*(*(Z_mol+0)-(*(Z_mol+1)));
		RMYW=sqrt(RMYW);
		RMYWA=RMYWA+RMYW;
		RTY=RMYWA/mm;
		cout<<mm<<"   "<<RTY<<endl;
	}*/


/*                                                      Enter into the main Loop                                                      */

	int NSTART=1;
	int NSTEP=0;
	while(NSTART<=MAXKB||KB<=MAXKB)
	{
		for(int NTIMES=NSTART;NTIMES<=MAXKB;NTIMES=NTIMES+1)
		{
			KB=KB+1;
			PREDCT(Tot_num);
			EVAL(Y,Nx,Ny,Nz,Tot_num,KB,KSORT,RDMAX,RDEL,RLIST,RMAX,TOTE,TOTV,NABTOT);
			CORR(Y,Nx,Ny,Nz,Tot_num,DELTA);
			TDIST=0;
			SUMVEL=0;
			for(int my=0;my=Tot_num;my=my+1)
			{
				TDIST=TDIST+(*(DISP_X+my))*(*(DISP_X+my))+(*(DISP_Y+my))*(*(DISP_Y+my))+(*(DISP_Z+my))*(*(DISP_Z+my));
				SUMVEL=SUMVEL+(*(XV_mol+my))*(*(XV_mol+my))+(*(YV_mol+my))*(*(YV_mol+my))+(*(ZV_mol+my))*(*(ZV_mol+my));
			}
	    	TDIST=TDIST/Tot_num;
	    	EK=SUMVEL/(2*Tot_num*DELTA*DELTA);

    		XSUM=XSUM+SUMVEL;
    		SUME=SUME+TOTE;
    		SUMV=SUMV+TOTV;

	    	if((KB%KWRITE)==0)
			{
	    		FKB=double(KB)*Tot_num;//All number which have been calculation.
	    		TMP=XSUM*(3*DELTA*DELTA*FKB);//temperature.
		    	ENR=(SUME/FKB)+CORE;//average potential energy on all have been calculation.
		    	VIR=VSCALE*(SUMV/FKB+CORV);
	     		PRES=TMP*DR*(1.0-VIR);
     			E1=TOTE/Tot_num+CORE;
      			ETOT=E1+EK;     //Total energy.
    			P1=TMP*DR*(1.0-VSCALE*(TOTV/Tot_num+CORV));
    			RLTIM=DELTA*KB*TSTEP;
     			if(KB%1000==0||KB==MAXKB)
    				RDF(DR,RDEL,NY);
    			if(IFLG<1)
					EQBRAT(KB,IFLG,AHEAT,SUMVEL,TDIST,XDIST,Tot_num,NY);
			}
		}
    	if(KB<=MAXKB)
    		NSTART=KB+1;
	}

/*                                                      Delete all arrays                                                      */

	delete []X_mol;
	delete []Y_mol;
	delete []Z_mol;

	delete []XV_mol;
	delete []YV_mol;
	delete []ZV_mol;	
	
	delete []XSE_mol;
	delete []YSE_mol;
	delete []ZSE_mol;
	
	delete []XTH_mol;
	delete []YTH_mol;
	delete []ZTH_mol;
	
	delete []XFO_mol;
	delete []YFO_mol;
	delete []ZFO_mol;

	delete []XFI_mol;
	delete []YFI_mol;
	delete []ZFI_mol;
	
	delete []ForX_mol;
	delete []ForY_mol;
	delete []ForZ_mol;

    delete []DISP_X;
    delete []DISP_Y;
    delete []DISP_Z;
	
	delete []XOLD_mol;
	delete []YOLD_mol;
	delete []ZOLD_mol;

	delete []NABORS;
	delete []IDIST;
    delete []LIST;
	
	return 1;
}