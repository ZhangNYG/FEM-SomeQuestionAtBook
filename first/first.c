//计算平面钢架的程序
//程序开始

#include<stdio.h>
#include<math.h>

#define NE 2
#define NJ 3                   //定义并输入基本参数
#define NZ 6
#define NPJ 3
#define NPF 2
#define NJ3 9
#define DD 9
#define E0 3.0000E7            // 定义并输入常数
#define A0 0.5
#define I0 4.16667E-2
#define PI 3.141592654

// 这是输入参数的初始化和定义全局变量

int jm[NE+1][3]={{0,0,0},{0,1,2},{0,3,1}};
double gc[NE+1]={0.0,5.0,5.0};
double gj[NE+1]={0.0,0.0,90.0};
double mj[NE+1]={0.0,A0,A0};
double gx[NE+1]={0.0,I0,I0};
int zc[NZ+1]={0,4,5,6,7,8,9};
double pj[NPJ+1][3]={{0.0,0.0,0.0},{0.0,6.0,1.0},{0.0,2.0,2.0},{0.0,-5.0,3.0}};
double pf[NPF+1][5]={{0,0,0,0,0},{0,-4.8,5.0,1.0,1.0},{0.0,-8.0,2.5,2.0,2.0}};
double kz[NJ3+1][NJ3+1],p[NJ3+1];
double pe[7],f[7],f0[7],t[7][7];
double ke[7][7],kd[7][7];
//  kz[][]---整体刚度矩阵
//  ke[][]---整体坐标系下的单元刚度矩阵
//  kd[][]---局部坐标系下的单元刚度矩阵
//  t[][] ---坐标板换矩阵



//*****函数声明*****
void jdugd(int);
void zb(int);
void gdnl(int);
void dugd(int);


//***主函数开始*****
void main(void)
{

	int i,j,k,e,dh,h,ii,jj,hz,a1,b1,m,l,dl,zl,z,j0;
	double cl,wy[7];
	int IM,IN,jn;


	//<功能：形成矩阵p>
	if(NPJ>0)
	{
		for(i=1;i<=NPJ;i++)
		{
			j=(int)pj[i][2];
			p[j]=pj[i][1];
		}
	}
	if(NPF>0)
	{
		for(i=1;i<=NPF;i++)
		{
			hz=i;
			gdnl(hz);
			e=(int)pf[hz][3];
			zb(e);
			for(j=1;j<=6;j++)
			{
				pe[j]=0.0;
				for(k=1;k<=6;k++)
				{
					pe[j]=pe[j]-t[k][j]*f0[k];    //用的是转秩矩阵     zb()函数出来的是直接矩阵   在这先让列不动行动  相当于乘的转秩
				}
			}
			a1=jm[e][1];
			b1=jm[e][2];
			p[3*a1-2]=p[3*a1-2]+pe[1];
			p[3*a1-1]=p[3*a1-1]+pe[2];
			p[3*a1]=p[3*a1]+pe[3];
			p[3*b1-2]=p[3*b1-2]+pe[4];
			p[3*b1-1]=p[3*b1-1]+pe[5];
			p[3*b1]=p[3*b1]+pe[6];

		}
	}






//*********************************************
	//<功能：生成整体刚度矩阵ke[][]>
	for(e=1;e<=NE;e++)
	{

		dugd(e);
		for(i=1;i<=2;i++)
		{
			for(ii=1;ii<=3;ii++)
			{
				h=3*(i-1)+ii;
				dh=3*(jm[e][i]-1)+ii;
				
			
		
				for(j=1;j<=2;j++)
				{
					for(jj=1;jj<=3;jj++)
					{
						l=3*(j-1)+jj;
						zl=3*(jm[e][j]-1)+jj;
						dl=zl-dh+1;
						if(dl>0)
							kz[dh][dl]=kz[dh][dl]+ke[h][l];
					}
				}
			}
		}
	}

	//****引入边界条件*******
	for(i=1;i<=NZ;i++)
	{
		z=zc[i];
		kz[z][1]=1.0;
		for(j=2;j<=DD;j++)
			kz[z][j]=0.0;
		if((z!=1))
		{
			if(z>DD)j0=DD;
			else if(z<=DD)j0=z;
			for(j=2;j<=j0;j++)
				kz[z-j+1][j]=0.0;
		}
		p[z]=0.0;                               //有点问题？？？？？？？？？？？？？！！！！！！！！
	}



	/*高斯消元法解方程组*/
	///



	//消元
	for(k=1;k<=NJ3-1;k++)
	{
		if(NJ3>k+DD-1)
			IM=k+DD-1;
		else if(NJ3<=k+DD-1)IM=NJ3;
		IN=k+1;                                   //you有问题
		for(i=IN;i<=IM;i++)
		{
			l=i-k+1;
			cl=kz[k][l]/kz[k][1];
			jn=DD-l+1;
			for(j=1;j<=jn;j++)
			{
				m=j+i-k;
				kz[i][j]=kz[i][j]-cl*kz[k][m];         //有问题

			}
			p[i]=p[i]-cl*p[k];

		}


	}

//******回代***********
	p[NJ3]=(p[NJ3]/kz[NJ3][1]);
	for(i=NJ3-1;i>=1;i--)
	{
		if(DD>NJ3-i+1)j0=NJ3-i+1;
		else j0=DD;
		for(j=2;j<=j0;j++)
		{
			h=j+i-1;
			p[i]=(p[i]-kz[i][j]*p[h]);
		}
		
		
			p[i]=(p[i]/kz[i][1]);
		

	}
	printf("\n");
	printf("\t\t有限元法（李景湧）第二章程序 \n");
	printf("__________________________________________________________\n");
	printf("NJ=		   U=	         V=	   	  CETA=	 \n");
	for(i=1;i<=NJ;i++)
	{
		printf("%-9d    %-12.11f   %-12.11f    %-12.11f\n",i,p[3*i-2],p[3*i-1],p[3*i]);
	}
	printf("__________________________________________________________\n");
	//*根据E的值输出相应E单元的N,Q,M（A,B）的结果**
	printf("E=   \t    N=        Q=           M=      \n");
	//计算轴力N‘,剪力Q，弯矩M*
		
	for(e=1;e<=NE;e++)
	{
		jdugd(e);
		zb(e);
		for(i=1;i<=2;i++)
		{
			for(ii=1;ii<=3;ii++)
			{
				h=3*(i-1)+ii;
				dh=3*(jm[e][i]-1)+ii;
				wy[h]=p[dh];


			}
		}
		for(i=1;i<=6;i++)
		
		{
			f[i]=0.0;
			for(j=1;j<=6;j++)
			{
				for(k=1;k<=6;k++)
				{
					f[i]=f[i]+kd[i][j]*t[j][k]*wy[k];
				}
			}
		}
		if(NPF>0)
		{
			for(i=1;i<=NPF;i++)
				if(pf[i][3]==e)
				{
					hz=i;
					gdnl(hz);
					for(j=1;j<=6;j++)
					{
						f[j]=f[j]+f0[j];
					}
				}
		}
		printf("%-4d(A)   %-9.5f  %-9.5f   %-9.5f\n",e,f[1],f[2],f[3]);
		printf("    (B)   %-9.5f  %-9.5f   %-9.5f \n",f[4],f[5],f[6]);


	}
	return;

}

//**************主程序结束**********************************



//gdnl()函数：<功能：将非节点荷载下的杆端力计算出来存入f0[]>
//***********************************************************

void gdnl(int hz)
{

	int ind,e;
	double g,c,l0,d;
	 g=pf[hz][1];
	 c=pf[hz][2];
	 e=(int)pf[hz][3];
	 ind=(int)pf[hz][4];
	 l0=gc[e];
	 d=l0-c;
	 if(ind==1)
	 
	 {
		 f0[1]=0.0;
		 f0[2]=-(g*c*(2-2*c*c/(l0*l0)+(c*c*c)/(l0*l0*l0)))/2;
		 f0[3]=-(g*c*c)*(6-8*c/l0+3*c*c/(l0*l0))/12;
		 f0[4]=0.0;
		 f0[5]=-g*c-f0[2];
		 f0[6]=(g*c*c*c)*(4-3*c/l0)/(12*l0);


	 }
	 else 
	 {
		 if(ind==2)
		 {
			 f0[1]=0.0;
			 f0[2]=(-(g*d*d)*(l0+2*c))/(l0*l0*l0);
			 f0[3]=-(g*c*d*d)/(l0*l0);
			 f0[4]=0.0;
			 f0[5]=(-(g*c*c)*(l0+2*d))/(l0*l0*l0);     //不知道书上是不是错了？
			 f0[6]=(g*c*c*d)/(l0*l0);

		 }
		 else
		 {
			 f0[1]=-(g*d/l0);
			 f0[2]=0.0;
			 f0[3]=0.0;
			 f0[4]=-g*c/l0;
			 f0[5]=0.0;
			 f0[6]=0.0;
		 }
	 }
}



//****************************************
//zb()函数：<功能：构成坐标变换矩阵>
//****************************************void
void zb(int e)
{
	double ceta,co,si;
	int i,j;
	ceta=(gj[e]*PI)/180;
	co=cos(ceta);
	si=sin(ceta);
	t[1][1]=co;
	t[1][2]=si;
	t[2][1]=-si;
	t[2][2]=co;
	t[3][3]=1.0;
	for(i=1;i<=3;i++)
	{
		for(j=1;j<=3;j++)
		{t[i+3][j+3]=t[i][j];}
	}

}







//****************************************************
//jdugd()函数：<功能：计算局部坐标下单元刚度矩阵kd[][]>
//*****************************************************
void jdugd(int e)
{
	double a0,l0,j0;
	int i,j;
	a0=mj[e];
	l0=gc[e];
	j0=gx[e];

	for(i=0;i<=6;i++)
		for(j=0;j<=6;j++)
			kd[i][j]=0.0;
	kd[1][1]=E0*a0/l0;
	kd[2][2]=12*E0*j0/pow(l0,3);
	kd[3][2]=6*E0*j0/pow(l0,2);
	kd[3][3]=4*E0*j0/l0;
	kd[4][1]=-kd[1][1];
	kd[4][4]=kd[1][1];
	kd[5][2]=-kd[2][2];
	kd[5][3]=-kd[3][2];
	kd[5][5]=kd[2][2];
	kd[6][2]=kd[3][2];
	kd[6][3]=2*E0*j0/l0;
	kd[6][5]=-kd[3][2];
	kd[6][6]=kd[3][3];

	for(i=1;i<=6;i++)
		for(j=1;j<=i;j++)
			kd[j][i]=kd[i][j];

}





//*****************************************
//dugd()函数<功能：计算整体坐标下单元刚度矩阵kc[][]>
//*******************************************
void dugd(int e)
{
	int i,j,k,m;
	jdugd(e);
	zb(e);
	for(i=1;i<=6;i++)
	{
		for(j=1;j<=6;j++)
		{
			ke[i][j]=0.0;
			for(k=1;k<=6;k++)
			{
				for(m=1;m<=6;m++)
				{
					ke[i][j]=ke[i][j]+t[k][i]*kd[k][m]*t[m][j];
				}
			}

			

		}
	}
}


//end