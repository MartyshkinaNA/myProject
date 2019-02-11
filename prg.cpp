#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <mpi.h>
#include <sys/time.h>
#include <cmath>
#pragma pack(1)
#define THR 3
// Структура f_info отвечает за информацию о графическом файле
struct f_info
{
	unsigned char signature[2];	// сигнатура типа файла BM
	unsigned int sizefile;		// размер файла
	unsigned int reserved;		// зарезервировано
	unsigned int addr_offset;	// смещение начала изображения
};
// Структура pic_info содержит информацию об изображении, которое содержиться в файле
struct pic_info
{
	  unsigned int Size;        	//  размер заголовка
	  unsigned int Width;       	//  ширина изображения в пикселях
	  unsigned int Height;      	//   высота изображения в пикселях
	  unsigned short int  Planes;   //  число плоскостей изображения
	  unsigned short int  BitCount; //  бит/пиксел 1,4,8,24
	  unsigned int Compression; 	//  тип сжатия
	  unsigned int SizeImage;   	//  размер сжатого изображения
	  unsigned int XPelsPerMeter;	// горизонтальное разрешение
	  unsigned int YPelsPerMeter;	// вертикальное разрешение пикс/м
	  unsigned int ClrUsed;     	//  количество используемых цветов
	  unsigned int ClrImportant;	//  число "важных" цветов
};
// Функция d вычисляет модуль вектора градиента
// a,b - ширина и высота обрабатываемого изображения
// i,j - координаты пикселя, в котором вычисляется модуль вектора градиента
// xm,ym - размер окрестности пиксяля
// *m1 - массив значений яркости пикселя исходного изображения
int d(int i, int j, unsigned int a, unsigned int b, double *m1, int xm, int ym)
{
	double gx=0; 
	double gy=0;
	// Вычисление частных производных: gx - по оси x; gy - по оси y
	for(int l=ym/2; l>0; l--)
		for(int k=-xm/2; k<=xm/2; k++)
		{
			gx+=m1[(i+k)*b+(j+l)]-m1[(i+k)*b+(j-l)]; 
			gy+=m1[(i+l)*b+(j+k)]-m1[(i-l)*b+(j+k)];	
		}
	double df=sqrt(gx*gx+gy*gy); // df - модуль вектора градиента
	return df;
}
// Функция u вычисляет направление (угол) вектора градиента
// a,b - ширина и высота обрабатываемого изображения
// i,j - координаты пикселя, в котором вычисляется угол вектора градиента
// xm,ym - размер окрестности пиксяля
// *m1 - массив значений яркости пикселя исходного изображения
double u(int i, int j, unsigned int a, unsigned int b, double *m1, int xm, int ym)
{
	double gx=0;
	double gy=0;
	// Вычисление частных производных: gx - по оси x; gy - по оси y
	for(int l=ym/2; l>0; l--)
		for(int k=-xm/2; k<=xm/2; k++)
		{
			gx+=m1[(i+k)*b+(j+l)]-m1[(i+k)*b+(j-l)];		
			gy+=m1[(i+l)*b+(j+k)]-m1[(i-l)*b+(j+k)];
		}
	double a1=atan2(gy,gx); // a1 - угол вектора градиента
	return a1;
}
// Функция outline выделяет контур объектов на изображении
// a,b - ширина и высота обрабатываемого изображения
// xm,ym - размер окрестности пиксяля
// *m1 - массив значений яркости пикселя исходного изображения
// *m2 - массив значений яркости пикселя обработанного изображения (изображения с контуром)
// P - маска для фильтрации (в данном случае маска Превитта)
void outline(unsigned int a, unsigned int b, double *m1, double *m2, int P[], int xm, int ym)
{	
	int E=25, df, ug; // Е - заданный неотрицательный порог
	double A=(15*3.14159)/180; // А - заданный неотрицательный угловой порог
	double sum;
	double *m;
	m=new double[a*b];
	for(unsigned int i=0;i<a;i++)
		for(unsigned int j=0;j<b;j++)
			m[i*b+j]=0;
	// фильтрация исходного изображения m1 с помощью масски P
	for(unsigned int i=xm/2;i<a-xm/2;i++)
		for(unsigned int j=ym/2;j<b-ym/2;j++)
		{
			sum=0;
			int ir=0;
			for(int k=-xm/2; k<=xm/2; k++)
				for(int l=-ym/2; l<=ym/2; l++)
					sum+=P[ir]*m1[(i+k)*b+(j+l)];
			m[i*b+j]=(sum>0)?sum:0; // m - результат фильтрации
		}
	// проверка условий, которые определяют является ли пиксель контурной точкой
	for(unsigned int i=xm-1;i<a-xm+1;i++)
		for(unsigned int j=ym-1;j<b-ym+1;j++)
		{	
			df=d(i,j,a,b,m,xm,ym);
			ug=u(i,j,a,b,m,xm,ym);
			for(int k=-xm/2; k<=xm/2; k++)
				for(int l=-ym/2; l<=ym/2; l++)
					if(((df-d(i+k,j+l,a,b,m,xm,ym))<=E)&&((ug-u(i+k,j+l,a,b,m,xm,ym))<=A))
						m2[(i+k)*b+(j+l)]=255; // m2 - финальный результат обработки
		}
}
// Функции time_start и time_stop ведут подсчет времени в миллисекундах
struct timeval tv1;
void time_start()
{
	struct timezone tz;
	gettimeofday(&tv1, &tz);
}
long time_stop()
{
	struct timeval tv2,dtv;
	struct timezone tz;
	gettimeofday(&tv2, &tz);
	dtv.tv_sec=tv2.tv_sec -tv1.tv_sec;
	dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
	if(dtv.tv_usec<0)
	{
		dtv.tv_sec--;
		dtv.tv_usec+=1000000;
	}
	return dtv.tv_sec*1000000+dtv.tv_usec;
}
int main(int argc, char *argv[])
{
	double *Bloc, *Bloc1;
	double *image, *ImOut;
	int P[9]={-1, -1, -1, 0, 0, 0, 1, 1, 1}; // маска Превитта
	int xm=3, ym=3;
	struct f_info f_i;
	struct pic_info pic_i;
	int size, rank, chunk, averag;
	char infile[100]="ish.bmp"; // "ish.bmp" - исходное изображение (формат изображения bmp 24 бит)
	char outfile[100]="obr.bmp"; // "obr.bmp" - обработанное изображение
	if(argc>1) strcpy(infile,argv[1]);
	if(argc>2) strcpy(outfile,argv[2]);
	std::ifstream in(infile,std::ios::in|std::ios::binary);
	if(!in) { std::cout<<"File not found...\n"; exit(1);}
	in.read((char*) &f_i, sizeof(f_info));
	in.read((char*) &pic_i, sizeof(pic_info));
	int ln_str=(f_i.sizefile-54)/pic_i.Height;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // определяет номер процесса 
	MPI_Comm_size(MPI_COMM_WORLD, &size); // определяет общее число процессов
	time_start();
	// инициализация массива исходного изображения image:
	if(rank==0)
	{
		image=new double[pic_i.Height*pic_i.Width];
		std::cout<<"Width - "<<pic_i.Width<<", Height - "<<pic_i.Height<<std::endl;
		for(unsigned int i=0;i<pic_i.Height;i++)
		{
			for(unsigned int j=0;j<pic_i.Width;j++)
			{
				char r, g, b;
				in.get(b); // b - значение интенсивности синего цвета
				in.get(g); // g - значение интенсивности зеленого цвета
				in.get(r); // r - значение интенсивности красного цвета
				image[i*pic_i.Width+j]=r*0.299+g*0.587+b*0.114; // преобразование в один сигнал яркрсти (переход к черно-белому изображению)
			}
			for(unsigned int k=0;k<(ln_str-pic_i.Width*3);k++)
			{
				char d;
				in.get(d);
			}
		}
		ImOut=new double[pic_i.Height*pic_i.Width];
		in.close();
		chunk=pic_i.Height/size; 
		averag=pic_i.Height%size;
	} 
	//  MPI_Bcast рассылает данные всем процессам
	MPI_Bcast(P,xm*ym,MPI_INT,0,MPI_COMM_WORLD); 
	MPI_Bcast(&xm,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ym,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&chunk,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&averag,1,MPI_INT,0,MPI_COMM_WORLD);
	int *sc, *ds, *str;
	sc = new int[size];
	ds= new int[size];
	str = new int[size];
	// цикл, с помощью которого определяется размер фрагмента изображения для каждого процесса (изображение разбиваеться на горизонтальные полосы):
	for(unsigned int i=0;i<size;i++)
	{
		int S=chunk;
		if(i<averag)
			S++;
		if(i==0 || i==size-1)
		{
			sc[i]=(S+2)*pic_i.Width;
			str[i] = S+2;
		}				
		else 
 		{
			sc[i]=(S+4)*pic_i.Width;
			str[i]=S+4;
		}	
		if(i==0)
			ds[i]=0;
		else
			ds[i]=(i*S-1)*pic_i.Width;
	}  
	Bloc=new double[str[rank]*pic_i.Width]; // массив для фрагмента исходного изображени (полосы)
	Bloc1=new double[str[rank]*pic_i.Width]; // массив для результата обработки фрагмента исходного изображения
	// MPI_Scatterv пересылает (распределяет) фрагметны (полосы) избражения по процессам
        MPI_Scatterv(image, sc, ds, MPI_DOUBLE, Bloc, sc[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  	outline(str[rank], pic_i.Width, Bloc, Bloc1, P, xm, ym); // вызов функции выделения контура на изображении
	if(rank!=0)
		for(unsigned int i=0;i<str[rank]-1;i++)
			for(unsigned int j=0;j<pic_i.Width;j++)
				Bloc1[i*pic_i.Width+j]=Bloc1[(i+1)*pic_i.Width+j];
	// цикл, вычисляющий размер обработанных фрагментов	
	for(unsigned int i=0;i<size;i++)
	{
		int S=chunk;
		if(i<averag)
			S++;
		sc[i]=(S)*pic_i.Width;
		str[i] = S;
		ds[i]=i*S*pic_i.Width;	
	} 
	// MPI_Gatherv собирает фрагменты изображения в единое
	MPI_Gatherv(Bloc1, sc[rank], MPI_DOUBLE, ImOut, sc, ds, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// запись результата в файл "obr.bmp":
	if(rank==0)
	{
		std::ofstream out;
		out.open(outfile,std::ios::out|std::ios::binary);
	 	out.write((char*) &f_i, sizeof(f_info));
		out.write((char*) &pic_i, sizeof(pic_info));
		for (unsigned int i=0;i<pic_i.Height;i++)
		{
			for (unsigned int j=0;j<pic_i.Width;j++)
			{
				out.put((char)ImOut[i*pic_i.Width+j]);
				out.put((char)ImOut[i*pic_i.Width+j]);
				out.put((char)ImOut[i*pic_i.Width+j]);
			}
			for(unsigned int k=0;k<(ln_str-pic_i.Width*3);k++) out.put(0);
		}
		std::cout<<"Press any key...\n";
		std::cout<<"Width - "<<pic_i.Width<<", Height - "<<pic_i.Height<<std::endl;
		delete [] ImOut;
		delete [] image;
		std::cout<<"Time: "<<time_stop()<<std::endl;
	}
	delete [] Bloc;
	delete [] Bloc1;
	delete [] ds;
	delete [] sc;
	delete [] str;

	MPI_Finalize();
		
	return 0;
}
