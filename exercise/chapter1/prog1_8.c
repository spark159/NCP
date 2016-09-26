#include <stdio.h>
#include <math.h>
#define SWAP(x,y,t) (t=x, x=y, y=t)
int count=0;

void perm(int [], int, int);
void main(void) {
  int list[10]={0,1,2,3,4,5,6,7,8,9};
  perm(list, 5, 0);
}

 
void perm(int list [], int n, int st ){
  if (st == n) {
    count +=1; printf ("#%3d ", count);
    int i;
    for (i=0; i<n; i++) {
      printf ("%d", list[i]);
    }
    printf ("\n");
  }

  else{
    int i,t;
    for (i=st; i<n; i++) {
      SWAP(list[st],list[i],t);
      perm(list, n, st+1);
    }    
  }
}

