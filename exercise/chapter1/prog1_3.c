#include <stdio.h>
#include <math.h>
#define MAX_SIZE 101
#define SWAP(x,y,t)(t=x,x=y,y=t)
void sort (int list[], int n);
void main (void)
{
  
  int i, n;
  printf("select length of list(should be less than %d):", MAX_SIZE);
  scanf("%d",&n);

  int list[n];
    
    for (i=0; i<n; i++) {
      list[i]=rand() % 1000;
    }

    sort(list, n);

    for (i=0; i<n; i++) {
      printf("%d", list[i]);
      printf("\n");
    }
}
  
void sort (int list[], int n)
{
  int i,j,min,t;
  for (i=0; i<n-1; i++) {
    min = list[i];
    for (j=i; j<n; j++) {
      if (list[j] < min){
	min=list[j]; SWAP(list[i],list[j],t);
      }
    }
  }
}
