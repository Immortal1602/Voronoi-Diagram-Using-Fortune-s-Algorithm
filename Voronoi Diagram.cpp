#include<bits/stdc++.h>
using namespace std;
#define n 1
#define Ep 0.000001

typedef struct point{
    double x;
    double y;
}point;


struct edge;


typedef struct face{
    point p;
    edge* e[n] = {NULL};
    int numOfEdges;

    face(point pp)
      : p(pp), numOfEdges(0) {}
    
}face;


typedef struct edge{
    point origin;
    point dest;
    bool done;
    face* f;
    edge * twin;
    // edge * next;
    // edge * prev;

    edge(point p)
      : origin(p), done(true) {}
    
}edge;


class Node{

    public:
    face* f;
    face* f1;
    face* f2;
    bool leaf;
    Node* left;
    Node* right;
    Node* parent;

    Node(face* ff){
        f = ff;
        f1 = NULL;
        f2 = NULL;
        leaf = true;
        left = NULL;
        right = NULL;
        parent = NULL;
    }
};


Node* root = NULL;


class event{

    public:
    point p;        // site point in case of site event and lowest most point of circle in case of circle event
    bool siteEvent; // site event or circle event
    Node* middleNode;
    Node* nextNode;
    Node* prevNode;
    bool valid;

    event()
    {
        siteEvent = true;
        middleNode = NULL;
        nextNode = NULL;
        prevNode = NULL;
        valid = true;
    }

    event(point pp, bool ss, Node* mn, Node* nn, Node* pn)
    {
        this -> p = pp;
        this -> siteEvent = ss;
        this -> middleNode = mn;
        this -> nextNode = nn;
        this -> prevNode = pn;
        this -> valid = true;
    }
};


class pq{

    public:
    event** arr;
    int capacity;
    int size;

    pq(int capacity)
    {
        arr = new event*[capacity];
        size = 0;
    }

    bool isEmpty(){
        return !(size);
    }

    void insert(point p, bool siteEvent, Node* mn, Node* nn, Node* pn){ // O(logn)
        size++;

        int index = size;
        arr[index] = new event(p, siteEvent, mn , nn, pn);

        while(index > 1)
        {
            int parent = index/2;

            if(arr[parent] -> p.y < arr[index] -> p.y)
            {
                swap(arr[parent], arr[index]);
                index = parent;
            }
            else
            {
                return;
            }

        }
    }

    void deletion(Node* node){  
        
        for(int i = 0 ; i < size ; i++){
            
            if(arr[i]->prevNode == node || arr[i]->middleNode == node || arr[i]->nextNode == node){

                arr[i]-> p.y = DBL_MAX;

                while(i != 1){
                    int parent = i/2;
                    swap(arr[parent] , arr[i]);
                    i = parent;
                }

                pop();
            }

        }
        
    }

    event* pop(){  // always root element is deleted   // O(logn)
        
        if(size == 0)
        {
            return NULL;
        }

        event* temp = arr[1];
        arr[1] = arr[size];
        size--;

        int i = 1;
        int largest = i;

        while(i<size)
        {
            int left = 2*i;
            int right = 2*i + 1;

            if(left <= size && arr[largest] -> p.y < arr[left] -> p.y)
            {
                largest = left;
            }
            if(right <= size && arr[largest] -> p.y < arr[right] -> p.y)
            {
                largest = right;
            }

            if(i != largest)
            {
                swap(arr[i],arr[largest]);
                i = largest;    
            }
            else
            {
                return temp;
            }
        }
    }
};

pq pq_events(3*n);


point intersection(point p0, point p1, double l)
{
    point res, p = p0;

    if (abs(p0.y - p1.y)<Ep)
        res.x = (p0.x + p1.x) / 2;
    else if (abs(p1.y - l)<Ep)
        res.x = p1.x;
    else if (abs(p0.y - l)<Ep) {
        res.x = p0.x;
        p = p1;
    } 
    else {
        // Use the quadratic formula.
        double z0 = 2*(p0.y - l);
        double z1 = 2*(p1.y - l);

        double a = 1/z0 - 1/z1;
        double b = -2*(p0.x/z0 - p1.x/z1);
        double c = ((p0.x*p0.x + p0.y*p0.y - l*l)/z0) - ((p1.x*p1.x + p1.y*p1.y - l*l)/z1);

        res.x = ( -b - sqrt(b*b - 4*a*c) ) / (2*a);
    }

    // Plug back into one of the parabola equations.
    res.y = (p.y*p.y + (p.x-res.x)*(p.x-res.x) - l*l)/(2*p.y-2*l);

    return res;
}

Node* insertBL(face* f){

    Node* temp = new Node(f);

    if(root == NULL)
    {
        root = temp;
        return NULL;
    }

    if(root -> leaf == 1){
        Node* temp2 = new Node(NULL); 

        if(root -> f -> p.x < temp -> f -> p.x){
            temp2 -> left = root;
            root = temp2;
            temp2 -> left -> parent = temp2;
            temp2 -> right = temp;
            temp -> parent = temp2;
            temp2 -> f1 -> p = root -> left -> f -> p;
            temp2 -> f2 -> p = temp -> f -> p;
            return temp2 -> left;
        }
        else{
            temp2 -> right = root;
            root = temp2;
            temp2 -> right -> parent = temp2;
            temp2 -> left = temp;
            temp -> parent = temp2;
            temp2 -> f1 -> p = temp -> f -> p;
            temp2 -> f2 -> p = root -> right -> f -> p;
            return temp2 -> right;
        }

    }

    Node* temp1 = root;

    while(temp1 -> leaf == 0)
    {
        point int_point = intersection(temp1 -> f1 -> p, temp1 -> f2 -> p, f -> p.y);
        if(temp -> f ->p.x < int_point.x)
            temp1 = temp1 -> left;
        else
            temp1 = temp1 -> right;
    }
    face* fj = temp1 -> f;

    Node* node1 = new Node(NULL);
    Node* node2 = new Node(NULL);
    Node* node3 = new Node(fj);
    Node* node4 = new Node(f);
    Node* node5 = new Node(fj);

    node1 -> leaf = node2 -> leaf = false;

    node1 -> f1 = fj;
    node1 -> f2 = f;

    node2 -> f1 = f;
    node2 -> f2 = fj;

    node3 -> parent = node1;

    node1 -> left = node3;
    node1 -> right = node2;
    node1 -> parent = temp1 -> parent;
    if(temp1 -> parent -> left == temp1)
        temp1 -> parent -> left = node1;
    else
        temp1 -> parent -> right = node1;

    node2 -> parent = node1;
    node2 -> right = node5;
    node2 -> left = node4;

    node4 -> parent = node2;

    node5 -> parent = node2;

    return node3;

}

void updateTuples(Node* node, face* f){

    if(node -> left -> leaf == 0 && node -> right -> leaf == 1){
        node -> f1 = node -> left -> f2;
        node -> f2 = node -> right -> f;
    }
    else if(node -> left -> leaf == 1 && node -> right -> leaf == 0){
        node -> f1 = node -> left -> f;
        node -> f2 = node -> right -> f1;
    }
    else if(node -> left -> leaf == 1 && node -> right -> leaf == 1){
        node -> f1 = node -> left -> f;
        node -> f2 = node -> right -> f;
    }
    else{
        node -> f1 = node -> left -> f2;
        node -> f2 = node -> right -> f1;
    }

    while(node -> parent != NULL){
        node = node -> parent;

        if((abs(node -> f1 -> p.x - f -> p.x)>Ep || (node -> f1 -> p.y - f -> p.y)>Ep) && (abs(node -> f2 -> p.x - f -> p.x )>Ep|| abs(node -> f2 -> p.y - f -> p.y)>Ep)){
            break;
        }

        node -> f1 = node -> left -> f2;
        node -> f2 = node -> right -> f1;
    }

}

void deleteBL(Node* delNode){
    Node* startUpdate = NULL;

    if(root -> leaf == true){
        root = NULL;
        return;
    }

    if(delNode -> parent == root){
        if(root -> left == delNode){
            root = root -> right;
            return;
        }
        if(root -> right == delNode){
            root = root -> left;
            return;
        }
    }

    if(delNode -> parent -> left == NULL || delNode -> parent -> right == NULL){

        Node* temp = delNode -> parent -> parent;
        startUpdate = temp -> parent;

        if(temp -> left == delNode -> parent){

            if(temp -> parent -> left == temp){
                temp -> parent -> left = temp -> right;
                temp -> right -> parent = temp -> parent;
            }
            else{
                temp -> parent -> right = temp -> right;
                temp -> right -> parent = temp -> parent;
            }

        }
        else{

            if(temp -> parent -> left == temp){
                temp -> parent -> left = temp -> left;
                temp -> left -> parent = temp -> parent;
            }
            else{
                temp -> parent -> right = temp -> left;
                temp -> left -> parent = temp -> parent;
            }

        }
    }
    else{

        Node* temp = delNode -> parent -> parent;
        startUpdate = temp;

        if(temp -> left == delNode -> parent){

            if(delNode -> parent -> left == delNode){
                temp -> left = delNode -> parent -> right;
                delNode -> parent -> right -> parent = temp;
            }
            else{
                temp -> left = delNode -> parent -> left;
                delNode -> parent -> left -> parent = temp;
            }
            
        }
        else{

            if(delNode -> parent -> left == delNode){
                temp -> right = delNode -> parent -> right;
                delNode -> parent -> right -> parent = temp;
            }
            else{
                temp -> right = delNode -> parent -> left;
                delNode -> parent -> left -> parent = temp;
            }

        }
    }

    updateTuples(startUpdate, delNode -> f);
    
}

void inorderfunc(Node* temp, face* f, face* f1, Node** mn, Node** pn){

    if(temp == NULL){
        return;
    }

    inorderfunc(temp -> left, f, f1, mn, pn);

    if(temp -> f1 == f1 && temp -> f2 == f){

        if(temp -> left -> leaf && !temp -> right -> leaf){
            (*pn) = temp -> left;
            (*mn) = temp -> right -> left;
            return;
        }
        else if(temp -> right -> leaf && !temp -> left -> leaf){
            (*mn) = temp -> right;
            (*pn) = temp -> left -> right;
            return;
        }
        else if(temp -> left -> leaf && temp -> right -> leaf){
            (*mn) = temp -> right;
            (*pn) = temp -> left;
            return;
        }

    }

    inorderfunc(temp -> right, f, f1, mn, pn);
}

void searchNode(face* f, face* f1, face* f2, Node** mn, Node** pn, Node** nn){

    inorderfunc(root, f , f1, mn, pn);

    inorderfunc(root, f2, f, nn, mn);

}

void HandleSiteEvent(face* f){

    Node* temp = insertBL(f);
    face* f1 = temp -> f;
    if(temp == NULL){
        return;
    }

    pq_events.deletion(temp);

    point new_start;
    new_start.y = ((temp-> f -> p.y)*(temp-> f -> p.y) + ((temp-> f -> p.x)- f -> p.x)*((temp-> f -> p.x)- f-> p.x) - (f -> p.y)*(f -> p.y))/(2*(temp-> f -> p.y)-2*(f -> p.y));
    new_start.x = f -> p.x;

    edge* e1 = new edge(new_start);
    edge* e2 = new edge(new_start);
    edge* e3 = new edge(new_start);
    edge* e4 = new edge(new_start);

    f->e[f->numOfEdges]=e1;
    f->e[f->numOfEdges]->twin=e2;
    f->numOfEdges += 1;

    f1->e[f1->numOfEdges]=e2;
    f1->e[f1->numOfEdges]->twin=e1;
    f1->numOfEdges += 1;


    f->e[f->numOfEdges]=e3;
    f->e[f->numOfEdges]->twin=e4;
    f->numOfEdges += 1;

    f1->e[f1->numOfEdges]=e4;
    f1->e[f1->numOfEdges]->twin=e3;
    f1->numOfEdges += 1;

    face* f2 = NULL;
    face* f3 = NULL;

    for(int i = 0; i < f1 -> numOfEdges; i++){

        if(!f1->e[i]->done){

            if(f1->e[i]->f->p.x < f1->p.x){

                f2 = f1 -> e[i] -> f;
                point temp = intersection(f1->p, f2->p, f->p.y);
                double dist1 = sqrt(pow(f1->e[i]->origin.x-new_start.x,2)+pow(f1->e[i]->origin.y-new_start.y,2));
                double dist2 = sqrt(pow(temp.y-new_start.y,2)+pow(temp.x-new_start.x,2));

                if(dist1 > dist2){
                    Node * middleNode = NULL , *nextNode = NULL , *prevNode = NULL; 
                    searchNode(f1,f2,f,&middleNode,&prevNode,&nextNode);

                    point a = f1 -> p;
                    point b = f2 -> p;
                    point c = f -> p;
                    double A = b.x - a.x,  B = b.y - a.y,
                    C = c.x - a.x,  D = c.y - a.y,
                    E = A*(a.x+b.x) + B*(a.y+b.y),
                    F = C*(a.x+c.x) + D*(a.y+c.y),
                    G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

                    point O;
                    O.x = (D*E-B*F)/G;
                    O.y = (A*F-C*E)/G;

                    point temp;
                    temp.x = O.x;
                    temp.y = O.x + sqrt( pow(a.x - O.x, 2) + pow(a.y - O.y, 2) );

                    pq_events.insert(temp , 0 , middleNode , nextNode , prevNode);
                }
            }
            else{

                f3 = f1 -> e[i] -> f;
                point temp = intersection(f1->p, f3->p, f->p.y);
                double dist1 = sqrt(pow(f1->e[i]->origin.x-new_start.x,2)+pow(f1->e[i]->origin.y-new_start.y,2));
                double dist2 = sqrt(pow(temp.y-new_start.y,2)+pow(temp.x-new_start.x,2));

                if(dist1 > dist2){
                    Node * middleNode = NULL , *nextNode = NULL , *prevNode = NULL; 
                    searchNode(f1,f,f3,&middleNode,&prevNode,&nextNode);
                    
                    point a = f1 -> p;
                    point b = f3 -> p;
                    point c = f -> p;
                    double A = b.x - a.x,  B = b.y - a.y,
                    C = c.x - a.x,  D = c.y - a.y,
                    E = A*(a.x+b.x) + B*(a.y+b.y),
                    F = C*(a.x+c.x) + D*(a.y+c.y),
                    G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

                    point O;
                    O.x = (D*E-B*F)/G;
                    O.y = (A*F-C*E)/G;

                    point temp;
                    temp.x = O.x;
                    temp.y = O.x + sqrt( pow(a.x - O.x, 2) + pow(a.y - O.y, 2) );

                    pq_events.insert( temp, 0 , middleNode , nextNode , prevNode);
                }
            }
        }
    }
}

void HandleCircleEvent(event * ev ,Node* node , face** faces){
    
    event* e = ev ;
    deleteBL(node);
    pq_events.deletion(node);

    point a = ev->middleNode->f->p;
    point b = ev->prevNode->f->p;
    point c = ev->nextNode->f->p;

   // Algorithm from O'Rourke 2ed p. 189.
    double A = b.x - a.x,  B = b.y - a.y,
           C = c.x - a.x,  D = c.y - a.y,
           E = A*(a.x+b.x) + B*(a.y+b.y),
           F = C*(a.x+c.x) + D*(a.y+c.y),
           G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

    point O;

    // Point o is the center of the circle.
    O.x = (D*E-B*F)/G;
    O.y = (A*F-C*E)/G;

    edge *e1 = new edge(O);
    edge *e2 = new edge(O);

    face* prevFace = e -> prevNode -> f1;
    face* nextFace = e -> nextNode -> f2;
    face * middleFace = e -> middleNode -> f;

    prevFace->e[prevFace->numOfEdges] = e1;
    prevFace->e[prevFace->numOfEdges]->f = prevFace;
    prevFace->e[prevFace->numOfEdges]->twin = e2;
    prevFace->numOfEdges++;

    nextFace->e[nextFace->numOfEdges] = e2;
    nextFace->e[nextFace->numOfEdges]->f = nextFace;
    nextFace->e[nextFace->numOfEdges]->twin = e1;
    nextFace->numOfEdges++;

    int i = -1;
    int j;
    while(middleFace->e[i] != NULL){
        i++;
        j = -1;

        while(nextFace->e[j] != NULL){
            j++;

            if(abs(nextFace->e[j]->origin.x - middleFace->e[i]->origin.x)<Ep && abs(nextFace->e[j]->origin.y == middleFace->e[i]->origin.y)<Ep){
                break;
            }

        }

    }

    nextFace->e[j]->dest = O;
    nextFace->e[j]->done = 1;
    middleFace->e[i]->dest = O;
    middleFace->e[i]->done = 1;

    i = -1;
    while(middleFace->e[i] != NULL){
        i++;
        j = -1;

        while(nextFace->e[j] != NULL){
            j++;

            if(abs(nextFace->e[j]->origin.x - middleFace->e[i]->origin.x)<Ep && abs(nextFace->e[j]->origin.y - middleFace->e[i]->origin.y)<Ep){
                break;
            }

        }

    }

    prevFace->e[j]->dest = O;
    prevFace->e[j]->done = 1;
    middleFace->e[i]->dest = O;
    middleFace->e[i]->done = 1;

    point mid = {(prevFace->p.x + nextFace->p.x)/2,(prevFace->p.y + nextFace->p.y)/2};

    for(int i = 0 ; i < prevFace->numOfEdges ; i++){

        if( !prevFace->e[i]->done ){

            point mid2 = {(prevFace->p.x + prevFace->e[i]->twin->f->p.x)/2,(prevFace->p.y + prevFace->e[i]->twin->f->p.y)/2};
            point origin2 = {prevFace->e[i]->origin.x , prevFace->e[i]->origin.y};

            if(sqrt(pow(mid.x-mid2.x,2) + pow(mid.y-mid2.y,2)) < sqrt(pow(O.x-origin2.x,2) + pow(O.y-origin2.y,2))){
                // NEED TO INSERT NODES INTO TREE AND THEN THE EVENT--THEN DO THE SAME THING FOR nextNode
                Node * middleNode = NULL , *nextNode = NULL , *prevNode = NULL; 
                searchNode(prevFace , prevFace->e[i]->f , nextFace  , &middleNode , &prevNode , &nextNode );

                point temp;
                temp.x = O.x;
                temp.y = O.x + sqrt( pow(a.x - O.x, 2) + pow(a.y - O.y, 2) );

                pq_events.insert(temp , 0 , middleNode , nextNode , prevNode);
            }

        }
    }

    for(int i = 0 ; i < prevFace->numOfEdges ; i++){

        if( !prevFace->e[i]->done ){

            point mid2 = {(prevFace->p.x + prevFace->e[i]->twin->f->p.x)/2,(prevFace->p.y + prevFace->e[i]->twin->f->p.y)/2};
            point origin2 = {prevFace->e[i]->origin.x , prevFace->e[i]->origin.y};

            if(sqrt(pow(mid.x-mid2.x,2) + pow(mid.y-mid2.y,2)) < sqrt(pow(O.x-origin2.x,2) + pow(O.y-origin2.y,2))){
                // NEED TO INSERT NODES INTO TREE AND THEN THE EVENT--THEN DO THE SAME THING FOR nextNode
                Node * middleNode = NULL , *nextNode = NULL , *prevNode = NULL; 
                searchNode(prevFace , prevFace->e[i]->f , nextFace  , &middleNode , &prevNode , &nextNode );

                point temp;
                temp.x = O.x;
                temp.y =  O.x + sqrt( pow(a.x - O.x, 2) + pow(a.y - O.y, 2) );

                pq_events.insert(temp , 0 , middleNode , nextNode , prevNode);
            }
        }
    }

 
}

double distance(point p1, point p2, int type) {

    switch(type) {

        case 1: return abs(p1.x - p2.x) + abs(p1.y - p2.y);
        
        case 2: return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));

        case 3: return max(abs(p1.x - p2.x), abs(p1.y - p2.y));

    }

}


int main(){

    point* pointSet = new point[n];

    face** faces = new face*[n];

    for(int i = 0 ; i < n ; i++){
        point p;
        cin >> p.x >> p.y;
        // pointSet[i] = p;
        faces[i] = new face(p);
        pq_events.insert(p, true, NULL, NULL, NULL);
    }

    while(!pq_events.isEmpty()){

        event* curr_ev = pq_events.pop();
        face* f;

        for(int i = 0 ; i < n; i++){
            if(abs((faces[i]-> p).x - (curr_ev -> p).x)<Ep && abs((faces[i]-> p).y - (curr_ev -> p).y)<Ep){
                f = faces[i];
                break;
            }
        }

        if(curr_ev -> siteEvent == 1){
            HandleSiteEvent(f);
        }

        else{
            HandleCircleEvent(curr_ev, curr_ev -> middleNode, faces);
        }
    }

    for(int i = 0 ; i <  n ; i ++){
        for(int j = 0 ; j < faces[i]->numOfEdges; j++){
            if(!faces[i]->e[j]->done){
                faces[i]->e[j]->done = 1;
                faces[i]->e[j]->dest = {(faces[i]->p.x + faces[i]->e[j]->twin->f->p.x)/2 , (faces[i]->p.y + faces[i]->e[j]->twin->f->p.y)/2};
            }
        }
    }

    cout<<"Enter x and y corrdinates of the point to be located : ";

    point pt;
    cin>>pt.x>>pt.y;


    int i;
    for( i = 0 ; i <  n ; i ++){
        int flag = 0;

        for(int j = 0 ; j < faces[i]->numOfEdges; j++){

            point focus = faces[i]->p;
            double m = (faces[i]->e[j]->origin.y -faces[i]->e[j]->dest.y)/(faces[i]->e[j]->dest.x -faces[i]->e[j]->origin.x) ;
            double c = (faces[i]->e[j]->origin.x*faces[i]->e[j]->dest.y - faces[i]->e[j]->dest.x*faces[i]->e[j]->origin.y)/(faces[i]->e[j]->dest.x -faces[i]->e[j]->origin.x) ;
            
            if((focus.y + m*focus.x + c)*(pt.y + m*pt.x + c) < 0){
                flag = 1;
                break;
            }

        }

        if(!flag){
            cout<<"The point lies in the voronoi cell of the point ("<<faces[i]->p.x<<","<<faces[i]->p.y<<")"<<endl;
            break;
        }

    }

    int nearest1 = 0;
    int nearest2 = 0;
    int nearest3 = 0;
    double minDist1 = DBL_MAX;
    double minDist2 = DBL_MAX;
    double minDist3 = DBL_MAX;

    for(int i = 0 ; i < n ; i++){
        
        double d1 = distance(faces[i] -> p, pt, 1);
        double d2 = distance(faces[i] -> p, pt, 2);
        double d3 = distance(faces[i] -> p, pt, 3);

        if(minDist1 > d1) {
            minDist1 = d1;
            nearest1 = i;
        }
        if(minDist2 > d2) {
            minDist2 = d2;
            nearest2 = i;
        }
        if(minDist3 > d3) {
            minDist3 = d3;
            nearest3 = i;
        }

    }

    cout<<"Nearest point according to Manhattan (L1) Distance is ("<<faces[nearest1]->p.x<<","<<faces[nearest1]->p.y<<")"<<endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    cout<<"Nearest point according to Euclidean (L2) Distance is ("<<faces[nearest2]->p.x<<","<<faces[nearest2]->p.y<<")"<<endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    cout<<"Nearest point according to Max Norm (L-infinity) Distance is ("<<faces[nearest3]->p.x<<","<<faces[nearest3]->p.y<<")"<<endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    
    return 0;
}