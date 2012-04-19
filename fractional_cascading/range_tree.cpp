#include <iostream>
#include <cstdlib>
#include <algorithm>
#define EPSILON 0.00000000001
#define LEFT_NODE 0
#define RIGHT_NODE 1
using namespace std;

struct cascade_node
{
	pair<double, double> p;
	struct cascade_node * left;
	struct cascade_node * right;
	int left_index;
	int right_index;
};

struct node
{
	pair<double, double> value;
	struct node* left;
	struct node* right;
	struct cascade_node *y_root;
	int Ny;
};

int debug =0;
void merge(struct cascade_node *A, int La, struct cascade_node * B, int Lb, struct cascade_node* merged);
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
			    double y1, double y2, int checkAncestor, int type, 
			    struct cascade_node* parent_low, struct cascade_node* parent_high);
int findRightSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    double y1, double y2, int checkAncestor, int type, 
			    struct cascade_node* parent_low, struct cascade_node* parent_high);
int inBoundedBox(double x, double y, double x1, double x2, double y1, double y2);
int binarySearch( double value, struct cascade_node a[], int l, int r );
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
	cout<<"Preorder \n";
	preorder(root);
*/
	for(i=0;i<R;i++)
	{
	//	cout<<"Range: "<<r_x[i].first<<" "<<r_x[i].second<<" "<<r_y[i].first<<" "<< r_y[i].second<<endl;
		
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
	
	struct cascade_node* y_subtree;
	
	int mid = (high + low)/2;
	struct node* root;
	root = new node;
	root->value = p_sx[mid];
	root->left =  preprocessTree(low, mid - 1);
	root->right = preprocessTree(mid + 1, high);
	
	if(root->left == NULL && root->right == NULL) // A leaf node
	{
		y_subtree = (struct cascade_node*) malloc(sizeof(struct cascade_node));
		y_subtree->p = p_sx[mid];
		y_subtree->left = NULL;
		y_subtree->right = NULL;
		y_subtree->left_index = -1;
		y_subtree->right_index = -1;
		//cout<<"Allocated "<< p_sx[mid].second <<" to leaf "<< p_sx[mid].first<<" "<<p_sx[mid].second;
		root->Ny = 1;
		root->y_root = y_subtree;
		//cout<<"Allocated "<< *root->y_root <<" to leaf "<< p_sx[mid].first<<" "<<p_sx[mid].second;
	}
	else if (root->left ==NULL)
	{
		int i= 0;
		int y_tree_len = root->right->Ny + 1; // Current node's information
		y_subtree = (struct cascade_node*) malloc( y_tree_len * sizeof(struct cascade_node) );
		
		for(i = 0; i<y_tree_len - 1; i++)
		{
			y_subtree[i].p =  root->right->y_root[i].p;
			y_subtree[i].left =  NULL;
			y_subtree[i].right = &root->right ->y_root[i];
			y_subtree[i].left_index = -1;
			y_subtree[i].right_index = i;
		}
		
		i=0;
		
		while(y_subtree[i].p.second < root->value.second && i < (y_tree_len -1) )i++;
		
		
		struct cascade_node temp = y_subtree[i]; 
		y_subtree[i].p = root->value;
		
		if(i>0)
		{
			y_subtree[i].left = y_subtree[i-1].left;
			y_subtree[i].right = y_subtree[i-1].right;
			y_subtree[i].left_index = y_subtree[i-1].left_index;
			y_subtree[i].right_index = y_subtree[i-1].right_index;
		}
		else
		{
			y_subtree[i].left = temp.left;
			y_subtree[i].right = temp.right;
			y_subtree[i].left_index = temp.left_index;
			y_subtree[i].right_index = temp.right_index;
		}
		
		i++;
		
		while(i < y_tree_len)
		{
			struct cascade_node k = y_subtree[i];
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
		y_subtree = (struct cascade_node*) malloc( y_tree_len * sizeof(struct cascade_node) );
		
		for(i = 0; i<y_tree_len - 1; i++)
		{
			y_subtree[i].p =  root->left->y_root[i].p;
			y_subtree[i].left = &root->left->y_root[i];
			y_subtree[i].right = NULL;
			y_subtree[i].left_index = i;
			y_subtree[i].right_index = -1;
		}
		
		i=0;
		
		while(y_subtree[i].p.second < root->value.second && i < (y_tree_len -1) )i++;
		
		struct cascade_node temp = y_subtree[i]; 
		y_subtree[i].p = root->value;
		
		if(i>0)
		{
			y_subtree[i].left = y_subtree[i-1].left;
			y_subtree[i].right = y_subtree[i-1].right;
			y_subtree[i].left_index = y_subtree[i-1].left_index;
			y_subtree[i].right_index = y_subtree[i-1].right_index;
		}
		else
		{
			y_subtree[i].left = temp.left;
			y_subtree[i].right = temp.right;
			y_subtree[i].left_index = temp.left_index;
			y_subtree[i].right_index = temp.right_index;
		}	
		
		i++;
		
		while(i < y_tree_len)
		{
			struct cascade_node k = y_subtree[i];
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
		y_subtree = (struct cascade_node*) malloc( y_tree_len * sizeof(struct cascade_node) );
		
		//merge 
		merge(root->left->y_root, root->left->Ny, root->right->y_root, root->right->Ny, y_subtree);
		//Insert root->value in the merged array
		while(y_subtree[i].p.second < root->value.second && i < (y_tree_len -1) )i++;
		
		struct cascade_node temp = y_subtree[i]; 
		y_subtree[i].p = root->value;
		
		if(i>0)
		{
			y_subtree[i].left = y_subtree[i-1].left;
			y_subtree[i].right = y_subtree[i-1].right;
			y_subtree[i].left_index = y_subtree[i-1].left_index;
			y_subtree[i].right_index = y_subtree[i-1].right_index;
		}
		else
		{
			y_subtree[i].left = temp.left;
			y_subtree[i].right = temp.right;
			y_subtree[i].left_index = temp.left_index;
			y_subtree[i].right_index = temp.right_index;
		}	
		 
		
		i++;
		
		while(i < y_tree_len)
		{
			struct cascade_node k = y_subtree[i];
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

void merge(struct cascade_node *A, int La, struct cascade_node * B, int Lb, struct cascade_node* merged)
{
	int i=0, j=0, k=0;
	while(i < La && j< Lb)
	{
		if(double_lt(A[i].p.second,B[j].p.second)) 
		{
			merged[k].p = A[i].p;
			merged[k].left = &A[i];
			merged[k].right = &B[j];
			merged[k].left_index = i;
			merged[k].right_index = j;
			k++;
			i++;
			
		}
		else 
		{
			merged[k].p = B[j].p;
			merged[k].left = &A[i];
			merged[k].right = &B[j];
			merged[k].left_index = i;
			merged[k].right_index = j;
			k++;
			j++;
		}
	}	
	while(i < La)
	{ 
		merged[k].p = A[i].p;
		merged[k].left = &A[i];
		merged[k].right = NULL;
		merged[k].left_index = i;
		merged[k].right_index = -1;
		i++;
		k++;
	}
	
	
	while(j < Lb) 
	{
		merged[k].p = B[j].p;
		merged[k].left = NULL;
		merged[k].right = &B[j];
		merged[k].left_index = -1;
		merged[k].right_index = j;
		j++;
		k++; 
	}
}

int searchTree(struct node* root, double x1, double x2, double y1, double y2)
{
	struct node* lca;
	//cout<<"Searching Left leaf";
	struct node* left_leaf = searchLeftLeaf(root, x1);
	//cout<<"Searching right leaf";
	struct node* right_leaf = searchRightLeaf(root, x2);
	//cout<<"End";
	if(left_leaf != right_leaf)
		lca = LCA(root, left_leaf, right_leaf);
	else lca = root; // There can be only one path!
	/*cout<<"Left Leaf "<<left_leaf->value.first<<" "<<left_leaf->value.second<<endl;
	cout<<"Right Leaf "<<right_leaf->value.first<<" "<<right_leaf->value.second<<endl;
	cout<<"LCA "<<lca->value.first<<" "<<lca->value.second<<endl;
	*/
	int N_max = lca->Ny;
	int lca_low_index = binarySearch(y1, lca->y_root, 0, lca->Ny -1);
	//struct cascade_node* lca_high = binarySearch (y2,  
	//cout<<"LCA L"<<lca_low_index;
	struct cascade_node* lca_low = &lca->y_root[lca_low_index];
	int left_result = 0;
	int right_result = 0;
	//cout<<lca_low_index<<" NM "<<N_max<<endl;
	if(lca_low && lca_low_index < N_max)
	{
	left_result= findLeftSubtree(lca->left, lca, x1, x2, y1, y2, 0, LEFT_NODE, lca_low, NULL);
	//cout<<"Left Nodes " << left_result<<endl;
	right_result = findRightSubtree(lca->right, lca, x1, x2, y1, y2, 0, RIGHT_NODE, lca_low, NULL);
	//cout<<"Right Nodes " << right_result<<endl;
	}
	
	if( inBoundedBox(lca->value.first, lca->value.second, x1, x2, y1, y2) )
	{
		//cout<<"Root";		
		return left_result + right_result + 1;
	}
	
	else return left_result + right_result;
	
}

int findLeftSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    double y1, double y2, int checkAncestor, int type, 
			    struct cascade_node* parent_low, struct cascade_node* parent_high)
{
	//cout<<"Enter"<<endl;
	//cout<<parent_low;
	int count = 0;
	double x, y;
	if(!root){return 0;}
	struct cascade_node* my_low = NULL;
	
	//cout<<"Before my_low"<<endl;
	if(type ==LEFT_NODE)
	{
		my_low = parent_low->left;
	}
	
	else{
		//cout<<"Null";
	//cout<<"Before else "<<parent_low;
		my_low = parent_low->right;
	}
	//cout<<"my low "<<my_low<<endl;
	
	if(checkAncestor == 1)
	{
		if(type == LEFT_NODE)
		{
			// Search in the parent's y subtree the range of y values
			if(parent->right)
			{
				//int start = binarySearch(y1, parent->right->y_root, 0, parent->right->Ny - 1);
				//cout<<" inside parent->right"<<endl;
				int start = -1;
				
				if(parent_low)
				start = parent_low->right_index; // magic!!

				//cout<<" Before while"<<endl;
				while(start < parent->right->Ny && start!=-1)
				{
					
					x = parent->right->y_root[start].p.first;
					y = parent->right->y_root[start].p.second;
						
					if(double_gt(y, y2))
					{ 
						break;
					}
					if( inBoundedBox(x, y, x1, x2, y1, y2) )
					{
//						cout<<x<<" "<<y<<endl;
						count ++;
					}
					start ++;
				}
				//cout<<" After While"<<endl;
				
			}
		}
		
		//cout<<" Before parent"<<endl;
		x = parent->value.first;
		y = parent->value.second;
		
		if( inBoundedBox(x, y, x1, x2, y1, y2) )
		{
//			cout<<x<<" "<<y<<endl;
			count ++;
		}
		
		
	}
	
	
	if(!root->left && !root->right)
	{
			x = root->value.first;
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
//				cout<<x<<" "<<y<<endl;
				count ++;
			}
			
			return count;
	}
	
	if(double_gt(x1, root->value.first))
	{
	//	cout<<"Before root->right && my_low"<<endl;
		if(root->right && my_low)
		{
	//		cout<<"Before root->right in && my_low"<<root->right<<" root "<<root<<" my_low "<<my_low<<endl;
			
			int result =  count + findLeftSubtree(root->right, root, x1, x2, y1, y2, 1, RIGHT_NODE, my_low, NULL);
	//		cout<<"After root->right && my_low"<<endl;
			return result;
			
		}
		else 
		{
			x = root->value.first; 
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{		
			//cout<<x<<" "<<y<<endl;
					count ++;
			}
			return count;    //if root->right isnt true, root->left has to be true and the leaf too! 
		
		}			  // This is because the tree is balanced
	}
	else //equal case
	{
		if(root->left)
		{
			if(my_low)
			{
		//	cout<<"Before my_low"<<endl;
			int result =  count + findLeftSubtree(root->left, root, x1, x2, y1, y2, 1, LEFT_NODE, my_low, NULL);
		//	cout<<"After my_low"<<endl;
			return result;
			}
			else
			{
				x = root->value.first; 
				y = root->value.second;
				if( inBoundedBox(x, y, x1, x2, y1, y2) )
				{
				//	cout<<x<<" "<<y<<endl;
					count ++;
				}
				return count;
			}
		}
		//	return searchLeftLeaf(root->left, low);
		else 
		{
		//We need to check both for itself and leaf before returning
		//	cout<<"Before final check"<<endl;
			x = root->value.first; 
			y = root->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
			//	cout<<x<<" "<<y<<endl;
				count ++;
			}
			
			
			x = root->right->value.first;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
			y = root->right->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
			//	cout<<x<<" "<<y<<endl;
				count ++;
			}
			//			cout<<"After final check"<<endl;
			return count;
		}
	}
	
	
	
}



int findRightSubtree(struct node *root, struct node* parent, double x1, double x2, 
			    double y1, double y2, int checkAncestor, int type, 
			    struct cascade_node* parent_low, struct cascade_node* parent_high)
{
	int count = 0;
	double x, y;
//cout<<"Init Parent "<< parent->value.first << " " << parent->value.second<<endl;
	if(!root){return 0;}
	
	struct cascade_node* my_low;
	
	if(type ==LEFT_NODE)
	{
		my_low = parent_low->left;
	}
	
	else
		my_low = parent_low->right;
	
	if(checkAncestor == 1)
	{
		if(type == RIGHT_NODE)
		{
			// Search in the parent's y subtree the range of y values
			//cout<<"Parent "<< parent->value.first << " " << parent->value.second<<endl;
			if(parent->left)
			{
				//int start = binarySearch(y1, parent->left->y_root, 0, parent->left->Ny - 1);
				int start = -1;
				
				if(parent_low)
				start = parent_low->left_index;
				
				while(start < parent->left->Ny && start!=-1)
				{
					//cout<<" start "<< start<<" ";	
					x = parent->left->y_root[start].p.first;
					y = parent->left->y_root[start].p.second;
					//cout<<"Outside box " << x<<" "<< y<< endl;
						
					if(double_gt(y, y2))
					{ 
						//cout<<"Broken " << x<<" "<< y<< endl;
						break;
					}
					if( inBoundedBox(x, y, x1, x2, y1, y2) )
					{
					//	cout<<x<<" "<<y<<endl;
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
		///	cout<<x<<" "<<y<<endl;
			count ++;
		}
		
		
	}
	
	
	if(!root->left && !root->right)
	{
			x = root->value.first;
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
		//		cout<<x<<" "<<y<<endl;
				count ++;
			}
			
			return count;
	}
	
	if(double_gt(x2, root->value.first) || double_equal(x2, root->value.first))
	{
		if(root->right)
		{
			//cout <<" I am here with count = " << count <<endl;
			if(my_low)
			{
			int result =  count + findRightSubtree(root->right, root, x1, x2, y1, y2, 1, RIGHT_NODE, my_low, NULL);
			//cout <<" I am here with result = " << result << endl;
			return result;
			}
			else
			{
				x = root->value.first; 
				y = root->value.second;
			
				if( inBoundedBox(x, y, x1, x2, y1, y2) )
				{
		//		cout<<x<<" "<<y<<endl;
					count ++;
				}
				return count;
			}
		}
		else 
		{
			
				x = root->value.first; 
				y = root->value.second;
			
				if( inBoundedBox(x, y, x1, x2, y1, y2) )
				{
		//		cout<<x<<" "<<y<<endl;
					count ++;
				}
			
			x = root->left->value.first;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
			y = root->left->value.second;
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
		//		cout<<x<<" "<<y<<endl;
				count ++;
			}
			
			return count;
		}   
	}
	else
	{
		if(root->left && my_low)
		{
					
			int result =  count + findRightSubtree(root->left, root, x1, x2, y1, y2, 1, LEFT_NODE, my_low, NULL);
			return result;
		}

		else 
		{
			x = root->value.first; 
			y = root->value.second;
			
			if( inBoundedBox(x, y, x1, x2, y1, y2) )
			{
		//			cout<<x<<" "<<y<<endl;
					count ++;
			}
			
			return count;
		}
	}
	
	
	
}


int binarySearch( double value, struct cascade_node a[], int l, int r )
{
    long low=l;
    long high = l>(r+1)?l:(r+1);

    while(low<high)
    {
        long mid = (low +high)/2;
        if ( double_lt(value,a[mid].p.second) || double_equal(value, a[mid].p.second) ) high =mid;
        else  low= mid+1;
    }
    return high;
}


struct node* searchLeftLeaf(struct node* root, double low)
{
//	cout<<root->value.first<<" "<<root->value.second<<endl;
	
	if(!root->left && !root->right)return root; 
	
	if(double_gt(low, root->value.first))
	{
		if(root->right)
			return searchLeftLeaf(root->right, low);
		else
		{
//		 cout<<root->left->value.first<<" "<<root->left->value.second<<endl;
	
		 return root->left;    //if root->right isnt true, root->left has to be true and the leaf too! 
		}			  // This is because the tree is balanced
	}
	else //equal case
	{
		if(root->left)
			return searchLeftLeaf(root->left, low);
		else 
		{
//		cout<<root->right->value.first<<" "<<root->right->value.second<<endl;
	
		return root->right;  //if root->left isnt true, root->RIGHT has to be true and the leaf too! 
		}	
	}			  // This is because the tree is balanced
}

struct node* searchRightLeaf(struct node* root, double high)
{
//	cout<<root->value.first<<" "<<root->value.second<<endl;
	
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
//		cout<<root->left->value.first<<" "<<root->left->value.second<<endl;
		return root->left;    //if root->right isnt true, root->left has to be true and the leaf too! 
		}			  // This is because the tree is balanced
	}
	else 
	{
		if(root->left)
			return searchRightLeaf(root->left, high);
		else 
		{
//			cout<<root->right->value.first<<" "<<root->right->value.second<<endl;
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
		cout<<root->y_root[i].p.first<< " "<<root->y_root[i].p.second<<" | ";
	
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
