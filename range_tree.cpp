#include <iostream>
#include <cstdlib>
#include <algorithm>
#define EPSILON 0.00000000001
#define LEFT_NODE 0
#define RIGHT_NODE 1
using namespace std;


struct node
{
	pair<double, double> value;
	struct node* left;
	struct node* right;
	pair<double, double> *y_root;
	int Ny;
};

int debug =0;
void merge(pair<double, double> *A, int La, pair<double, double>* B, int Lb, pair<double, double>* merged);
bool sortX (pair<double,double> i,pair<double,double> j);
bool sortY (pair<double,double> i,pair<double,double> j);
int searchTree(struct node* root, double x1, double x2, double y1, double y2);
int double_equal(double a, double b);
int double_lt(double a, double b);
int double_gt(double a, double b);
void inorder(struct node* root);
void preorder(struct node* root);
struct node* preprocessTree(int low, int high);
struct node* searchLeftLeaf(struct node* root, double low);
struct node* searchRightLeaf(struct node* root, double high);
struct node *LCA(struct node *root, struct node *p, struct node *q);
int findLeftSubtree(struct node *root, struct node* parent, double x1, double x2,
			          double y1, double y2, int checkAncestor, int type);
int findRightSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    	  double y1, double y2, int checkAncestor, int type);
int inBoundedBox(double x, double y, double x1, double x2, double y1, double y2);
int binarySearch( double value, pair<double, double> a[], int l, int r );
pair<double,double> *p_sx; 
int main()
{
	//Takes in N points <x,y> and R <x1,x2,y1,y2> Ranges
	//x1<x2; y1<y2
	int N, R, i; 
	cin>>N>>R;

	p_sx = (pair<double,double>*)malloc(sizeof(pair<double,double>)*N);
	
	pair<double,double> r_x[R];
	pair<double,double> r_y[R];
	double a, b, c, d;
	//Input pointset
	for(i=0;i<N;i++)
	{
		cin>>p_sx[i].first>>p_sx[i].second;

	
	}
	
	for(i=0;i<R;i++)
	{
		cin>> a >> b >> c >> d;
		
		//Ensuring x1 < x2 & y1 < y2 
		// Handle degenrate case of x1 = x2 / y1 = y2?
		
		if(double_lt(a,b))
		{
			r_x[i].first = a;
			r_x[i].second = b;
		} 
		else
		{
			r_x[i].first = b;
			r_x[i].second = a;
		}
		
		if(double_lt(c,d))
		{
			r_y[i].first = c;
			r_y[i].second = d;
		} 
		else
		{
			r_y[i].first = d;
			r_y[i].second = c;
		}
		
	}
	//points sorted by the x coordinate
	
	sort(p_sx, p_sx+N, sortX);
	
	struct node* root;
	root = preprocessTree(0, N-1);
	free(p_sx);
	
/*	cout<<"Root\n"<<root->value.first<< " "<<root->value.second;
	cout<<"Inorder \n";
	inorder(root);
/*	cout<<"Preorder \n";
	preorder(root);
*/
	for(i=0;i<R;i++)
	{
	//	cout<<"Range: "<<r_x[i].first<<" "<<r_x[i].second<<" "<<r_y[i].first<<" "<< r_y[i].second<<endl;
		if(double_equal(r_x[i].first, 1.0) && double_equal(r_x[i].second, 8.0) && 
		double_equal(r_y[i].first, 37.0) && double_equal(r_x[i].first, 39.0) )
		debug =0;
		
		int num_points = searchTree(root, r_x[i].first, r_x[i].second, r_y[i].first, r_y[i].second);
		//cout<<"Total points = "<<num_points<<endl;
		cout<<num_points<<endl;
	}

}

struct node* preprocessTree(int low, int high)
{
	if(high < low)
	{
		return NULL;
	}
	
	pair<double, double>* y_subtree;
	
	int mid = (high + low)/2;
	struct node* root;
	root = new node;
	root->value = p_sx[mid];
	root->left =  preprocessTree(low, mid - 1);
	root->right = preprocessTree(mid + 1, high);
	
	if(root->left == NULL && root->right == NULL) // A leaf node
	{
		y_subtree = (pair<double, double>*) malloc(sizeof(double));
		*y_subtree = p_sx[mid];
		//cout<<"Allocated "<< p_sx[mid].second <<" to leaf "<< p_sx[mid].first<<" "<<p_sx[mid].second;
		root->Ny = 1;
		root->y_root = y_subtree;
		//cout<<"Allocated "<< *root->y_root <<" to leaf "<< p_sx[mid].first<<" "<<p_sx[mid].second;
	}
	else if (root->left ==NULL)
	{
		int i= 0;
		int y_tree_len = root->right->Ny + 1; // Current node's information
		y_subtree = (pair<double, double>*) malloc( y_tree_len * sizeof(pair<double, double>) );
		
		for(i = 0; i<y_tree_len - 1; i++)
		{
			y_subtree[i] =  root->right->y_root[i];
		}
		
		i=0;
		
		while(y_subtree[i].second < root->value.second && i < (y_tree_len -1) )i++;
		
		pair<double, double> temp = y_subtree[i]; 
		
		y_subtree[i] = root->value;
		i++;
		
		while(i < y_tree_len)
		{
			pair<double, double> k = y_subtree[i];
			y_subtree[i] = temp;
			temp = k;
			i++;
		}
		
		//y_root now points to the merged subtree
		root->y_root = y_subtree;
		root->Ny = y_tree_len;
	}
	else if (root->right ==NULL)
	{
		int i= 0;
		int y_tree_len = root->left->Ny + 1; // Current node's information
		y_subtree = (pair<double, double>*) malloc( y_tree_len * sizeof(pair<double, double>) );
		
		for(i = 0; i<y_tree_len - 1; i++)
		{
			y_subtree[i] =  root->left->y_root[i];
		}
		
		i=0;
		
		while(y_subtree[i].second < root->value.second && i < (y_tree_len -1) )i++;
		
		pair<double, double> temp = y_subtree[i]; 
		y_subtree[i] = root->value;
		i++;
		
		while(i < y_tree_len)
		{
			pair<double, double> k = y_subtree[i];
			y_subtree[i] = temp;
			temp = k;
			i++;
		}
		
		//y_root now points to the merged subtree
		root->y_root = y_subtree;
		root->Ny = y_tree_len;
	}
	else if(root->left && root->right)
	{
		int i= 0;
		int y_tree_len = root->left->Ny + root->right->Ny + 1; // Current node's information
		y_subtree = (pair<double, double>*) malloc( y_tree_len * sizeof(pair<double, double>) );
		
		//merge 
		merge(root->left->y_root, root->left->Ny, root->right->y_root, root->right->Ny, y_subtree);
		//Insert root->value in the merged array
		while(y_subtree[i].second < root->value.second && i < (y_tree_len -1) )i++;
		
		pair<double, double> temp = y_subtree[i]; 
		y_subtree[i] = root->value;
		i++;
		
		while(i < y_tree_len)
		{
			pair<double, double> k = y_subtree[i];
			y_subtree[i] = temp;
			temp = k;
			i++;
		}
		
		//y_root now points to the merged subtree
		root->y_root = y_subtree;
		root->Ny = y_tree_len;
	}
	
	return root;

}

//merges A and B into merged

void merge(pair<double, double> *A, int La, pair<double, double>* B, int Lb, pair<double, double>* merged)
{
	int i=0, j=0, k=0;
	while(i < La && j< Lb)
	{
		if(double_lt(A[i].second,B[j].second)) merged[k++] = A[i++];
		else merged[k++] = B[j++];
	}	
	while(i < La) merged[k++] = A[i++];
	while(j < Lb) merged[k++] = B[j++]; 
	
}

int searchTree(struct node* root, double x1, double x2, double y1, double y2)
{
	struct node* left_leaf = searchLeftLeaf(root, x1);
	//cout<<"Searching right leaf";
	struct node* right_leaf = searchRightLeaf(root, x2);
	//cout<<"End";
	struct node* lca = LCA(root, left_leaf, right_leaf);
	//cout<<"Left Leaf "<<left_leaf->value.first<<" "<<left_leaf->value.second<<endl;
	//cout<<"Right Leaf "<<right_leaf->value.first<<" "<<right_leaf->value.second<<endl;
	//cout<<"LCA "<<lca->value.first<<" "<<lca->value.second<<endl;
	
	int left_result = findLeftSubtree(lca->left, lca, x1, x2, y1, y2, 0, LEFT_NODE);
	//cout<<"Left Nodes " << left_result<<endl;
	int right_result = findRightSubtree(lca->right, lca, x1, x2, y1, y2, 0, RIGHT_NODE);
	//cout<<"Right Nodes " << right_result<<endl;
	if( inBoundedBox(lca->value.first, lca->value.second, x1, x2, y1, y2) )
	{
		
		return left_result + right_result + 1;
	}
	
	else return left_result + right_result;
	
}

int findLeftSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    double y1, double y2, int checkAncestor, int type)
{
	int count = 0;
	double x, y;
	if(!root){return 0;}
	if(checkAncestor == 1)
	{
		if(type == LEFT_NODE)
		{
			// Search in the parent's y subtree the range of y values
			if(parent->right)
			{
				int start = binarySearch(y1, parent->right->y_root, 0, parent->right->Ny - 1);
				if(debug)
				cout<<"Parent "<< parent->value.first << " " << parent->value.second<<endl;
				while(start < parent->right->Ny)
				{
					if(debug)
					cout<<" start "<< start<<" ";	
					
					x = parent->right->y_root[start].first;
					y = parent->right->y_root[start].second;
					if(debug)
					cout<<"Outside box " << x<<" "<< y<< endl;
						
					if(double_gt(y, y2))
					{ 
						if(debug)
						cout<<"Broken " << x<<" "<< y<< endl;
						break;
					}
					if( inBoundedBox(x, y, x1, x2, y1, y2) )
					{
						if(debug)
						cout<<"In box " << x<<" "<< y<< endl;
						count ++;
					}
					start ++;
				}
				
			}
		}
		
		x = parent->value.first;
		y = parent->value.second;
		
		if( inBoundedBox(x, y, x1, x2, y1, y2) )
		{
			count ++;
		}
		
		
	}
	
	
	if(!root->left && !root->right)
	{
			x = root->value.first;
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
				count ++;
			}
			
			return count;
	}
	
	if(double_gt(x1, root->value.first))
	{
		if(root->right)
		{
			
			int result =  count + findLeftSubtree(root->right, root, x1, x2, y1, y2, 1, RIGHT_NODE);
			return result;
		}
		else 
		{
			x = root->value.first; 
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
					count ++;
			}
			return count;    //if root->right isnt true, root->left has to be true and the leaf too! 
		
		}			  // This is because the tree is balanced
	}
	else //equal case
	{
		if(root->left)
		{
			int result =  count + findLeftSubtree(root->left, root, x1, x2, y1, y2, 1, LEFT_NODE);
			return result;
		}
		//	return searchLeftLeaf(root->left, low);
		else 
		{
		//We need to check both for itself and leaf before returning
			
			x = root->value.first; 
			y = root->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
				count ++;
			}
			
			
			x = root->right->value.first;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
			y = root->right->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
				count ++;
			}
			
			return count;
		}
	}
	
	
	
}



int findRightSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    double y1, double y2, int checkAncestor, int type)
{
	int count = 0;
	double x, y;
//cout<<"Init Parent "<< parent->value.first << " " << parent->value.second<<endl;
	if(!root){return 0;}
	
	if(checkAncestor == 1)
	{
		if(type == RIGHT_NODE)
		{
			// Search in the parent's y subtree the range of y values
			//cout<<"Parent "<< parent->value.first << " " << parent->value.second<<endl;
			if(parent->left)
			{
				int start = binarySearch(y1, parent->left->y_root, 0, parent->left->Ny - 1);

				while(start < parent->left->Ny)
				{
					//cout<<" start "<< start<<" ";	
					x = parent->left->y_root[start].first;
					y = parent->left->y_root[start].second;
					//cout<<"Outside box " << x<<" "<< y<< endl;
						
					if(double_gt(y, y2))
					{ 
						//cout<<"Broken " << x<<" "<< y<< endl;
						break;
					}
					if( inBoundedBox(x, y, x1, x2, y1, y2) )
					{
						//cout<<"In box " << x<<" "<< y<< endl;
						count ++;
					}
					start ++;
				}
				
			}
			
		}
		
		// We always need to check the parent
		x = parent->value.first;
		y = parent->value.second;
			
		if( inBoundedBox(x, y, x1, x2, y1, y2) )
		{
			count ++;
		}
		
		
	}
	
	
	if(!root->left && !root->right)
	{
			x = root->value.first;
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
				count ++;
			}
			
			return count;
	}
	
	if(double_gt(x2, root->value.first) || double_equal(x2, root->value.first))
	{
		if(root->right)
		{
			//cout <<" I am here with count = " << count <<endl;
			int result =  count + findRightSubtree(root->right, root, x1, x2, y1, y2, 1, RIGHT_NODE);
			//cout <<" I am here with result = " << result << endl;
			return result;
		}
		else 
		{
			
				x = root->value.first; 
				y = root->value.second;
			
				if( inBoundedBox(x, y, x1, x2, y1, y2) )
				{
					count ++;
				}
			
			x = root->left->value.first;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
			y = root->left->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
				count ++;
			}
			
			return count;
		}   
	}
	else
	{
		if(root->left)
		{
					
			int result =  count + findRightSubtree(root->left, root, x1, x2, y1, y2, 1, LEFT_NODE);
			return result;
		}

		else 
		{
			x = root->value.first; 
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
					count ++;
			}
			
			return count;
		}
	}
	
	
	
}


int binarySearch( double value, pair<double, double> a[], int l, int r )
{
    long low=l;
    long high = l>(r+1)?l:(r+1);

    while(low<high)
    {
        long mid = (low +high)/2;
        if ( double_lt(value,a[mid].second) || double_equal(value, a[mid].second) ) high =mid;
        else  low= mid+1;
    }
    return high;
}


struct node* searchLeftLeaf(struct node* root, double low)
{
	if(!root->left && !root->right)return root; 
	
	if(double_gt(low, root->value.first))
	{
		if(root->right)
			return searchLeftLeaf(root->right, low);
		else return root->left;    //if root->right isnt true, root->left has to be true and the leaf too! 
					  // This is because the tree is balanced
	}
	else //equal case
	{
		if(root->left)
			return searchLeftLeaf(root->left, low);
		else return root->right;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
	}				  // This is because the tree is balanced
}

struct node* searchRightLeaf(struct node* root, double high)
{
	//cout<<root->value.first<<" "<<root->value.second<<endl;
	
	if(!root->left && !root->right)
	{
		return root; 
	}
	if(double_gt(high, root->value.first) || double_equal(high, root->value.first))
	{
		if(root->right)
			return searchRightLeaf(root->right, high);
		else 
		{
		//cout<<root->left->value.first<<" "<<root->left->value.second<<endl;
		return root->left;    //if root->right isnt true, root->left has to be true and the leaf too! 
		}			  // This is because the tree is balanced
	}
	else 
	{
		if(root->left)
			return searchRightLeaf(root->left, high);
		else 
		{
		//	cout<<root->right->value.first<<" "<<root->right->value.second<<endl;
		return root->right;
		}  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
	}	
}


struct node *LCA(struct node *root, struct node *p, struct node *q) {
  if (!root || !p || !q) return NULL;
  if (max(p->value.first, q->value.first) < root->value.first)
    return LCA(root->left, p, q);
  else if (min(p->value.first, q->value.first)  > root->value.first)
    return LCA(root->right, p, q);
  else
    return root;
}

void inorder(struct node* root)
{
	if (!root) return;
	inorder(root->left);
	
	cout<< "X "<<root->value.first<<" Y "<<root->value.second<<endl;
	
	inorder(root->right); 
}

void preorder(struct node* root)
{
	if (!root) return;
	cout<< "X "<<root->value.first<<" Y "<<root->value.second<<endl;
	cout<<"Sorted Y Array \n";
	
	int len = root->Ny;
	for(int i=0 ; i<len ; i++)
		cout<<root->y_root[i].first<< " "<<root->y_root[i].second<<" | ";
	
	cout<<endl;
	preorder(root->left);	
	preorder(root->right); 
}

int double_equal(double a, double b)
{
	if( ((a-b) < EPSILON && (a-b) >= 0) || ((b-a) < EPSILON && (b-a) >= 0) )
		return 1;
	return 0;
	     
}

int double_lt(double a, double b)
{
	if((b - a) > EPSILON) return 1;
	return 0;
}


int double_gt(double a, double b)
{
	if((a - b) > EPSILON) return 1;
	return 0;
}

int inBoundedBox(double x, double y, double x1, double x2, double y1, double y2)
{
	if( (double_equal(x,x1) || double_gt(x, x1) ) && 
	    (double_equal(x,x2) || double_lt(x, x2) ) &&
	    (double_equal(y,y1) || double_gt(y, y1) ) && 
	    (double_equal(y,y2) || double_lt(y, y2) ) )
	    
	    return 1;
	    
	return 0;
}

bool sortX (pair<double,double> i,pair<double,double> j) 
{ 
	if(i.first < j.first)
		return true;
	else if (i.first == j.first && i.second < j.second)
		return true;
	else return false;	
}

bool sortY (pair<double,double> i,pair<double,double> j) 
{ 
	if(i.second < j.second)
		return true;
	else if (i.second == j.second && i.first < j.first)
		return true;
	else return false;	
}
