#include <iostream>
#include <algorithm>
using namespace std;
struct Node{
    double x, y;
    int ind;
}*Point, *CH;
int Monotone_Chain(int n);
double cross(Node o, Node a, Node b)   {return (a.x-o.x)*(b.y-o.y)-(b.x-o.x)*(a.y-o.y);}
bool cmp(Node a, Node b) {return a.x < b.x || (a.x == b.x && a.y < b.y);} //first cmp x, next cmp y
int main()
{
    int n;
    while(cin >> n){ //輸入方式要改成直接讀array
        Point = new Node[n];
        CH = new Node[n*2];
        for(int i = 0; i < n; i++){
            cin >> Point[i].x >> Point[i].y; //這裡要改成輸入x = mu, y = reliability
	        Point[i].ind = i+1;
	}
        int m = Monotone_Chain(n);
        for(int i = 0; i <= m; i++){     //輸出圍成凸包的點
            cout << CH[i].x << " " << CH[i].y << endl;
	    cout << CH[i].ind<< endl;
}
        delete [] Point;
        delete [] CH;
    }
    return 0;
}
int Monotone_Chain(int n){
    sort(Point, Point+n, cmp);
    int m = 0;

    for(int i = 0; i < n; i++){
        while(m >= 2 && cross(CH[m-2], CH[m-1], Point[i]) >= 0)   m--;
        CH[m++] = Point[i];
    }
  

    for(int i = n-2, t = m+1; i >= 0; i--){     //t = m+1, 因為已經包括終點了
        while(m >= t && cross(CH[m-2], CH[m-1], Point[i]) >= 0)   m--;
        CH[m++] = Point[i];
    }

    return m-1; //m-1去除起點, 所以總共有m-1個頂點
} 
