#include "ExpressionParserData.h"

using namespace std;


node* createVariableNode(char* nodeName){
	node* n;
	n = (node*) malloc( sizeof(node) );
	//set type as input
	n->type = 0; //TODO put flopoco namespace flopoco::FPPipeline::FPNode::input
	n->name = nodeName;
	return n;
}

node* createConstantNode(int nodeValue){
	node* n;
	n = (node*) malloc( sizeof(node) );
	//set type as input
	n->type = 0; //TODO put flopoco namespace flopoco::FPPipeline::FPNode::input
	n->value = nodeValue;
	return n;
}

node* createNodeHavingRightOperandOnly(node* rightNode, int opType){
	node* n;
	n = (node*) malloc( sizeof(node) );
	//set type as input
	n->type = opType;
	nodeList* t;
	t=(nodeList*)malloc( sizeof(nodeList));
	t->n = rightNode;
	t->next = NULL;
	n->nodeArray = t; 
	return n;
}

//operations on lists

nodeList* appendNodeToListBack(node* n, nodeList* nl){
	//nl is the list head, so store it before iterating using it
	nodeList* head = nl;
	//get to the last element from the list
	while (nl->next!=NULL)
		nl = nl->next;	
	//create the new list element we want to insert
	nodeList *tail;
	tail = (nodeList*) malloc( sizeof( nodeList) );
	tail->n = n; //copy the pointer
	tail->next = NULL; //we point to heaven
	nl->next = tail; // linking done
	//return the head of the list
	return head;
}

nodeList* appendNodeToListFront(node* n, nodeList* nl){
	//nl is the list head, so store it before iterating using it

	//create the new list element we want to insert
	nodeList *elem;
	elem = (nodeList*) malloc( sizeof( nodeList) );
	elem->n = n; //copy the pointer
	elem->next = nl; //we point to heaven
	//return the head of the list
	return elem;
}

nodeList* createNodeList(node* n){
	nodeList* head;
	head=(nodeList*) malloc( sizeof(nodeList) );
	head->n = n; //copy pointer value
	head->next = NULL;	
	return head;
}

varList* appendVarToList(char* s, varList* vl){
	//nl is the list head, so store it before iterating using it
	varList* head = vl;
	//get to the last element from the list
	while (vl->next!=NULL)
		vl = vl->next;	
	//create the new list element we want to insert
	varList *tail;
	tail = (varList*) malloc( sizeof(varList));
	tail->name = s; //copy the pointer
	tail->next = NULL; //we point to heaven
	vl->next = tail; // linking done
	//return the head of the list
	return head;
}

varList* createVarList(char* s){
	varList* head;
	head=(varList*) malloc( sizeof(varList));
	head->name = s; //copy pointer value
	head->next = NULL;	
	return head;
}

void printExpression(node* n){
	if (n->name!=NULL)
		cout << n->name << "=";
	cout <<  "op"<<n->type;
	cout << "[";
	nodeList* head = n->nodeArray;
	while (head!=NULL){
		if (head->n->type==0){ //input node
			if (head->n->name!=NULL) //variable
				cout << head->n->name;
			else //constant
				cout << head->n->value;	
		} else {
			cout << "(";
			printExpression(head->n);	
			cout << ")";
		}
		head = head->next;
	}
	cout << "]";
}

void makeComputationalTree(node* parent, nodeList* expressionList, nodeList* statementList){
	nodeList* h = expressionList;
	while (h!=NULL){ //all statements here
		if (h->n->type==0){ // if it's an input node
			//if this node is not a constant
			if (h->n->name!=NULL){
				//fetch the name;
				char* c = h->n->name;
//				cout << "Found input node with name:" << c << endl;
				//parse	statements to see if const is declared
				nodeList* th = statementList;
				bool found=false;
				while ((th!=NULL)&&(!found)){
					if ( strcmp(th->n->name,c) ==  0)
						found = true;	
					else
					th = th->next;	
				}
				if (th!=NULL){ //found declaration: replace constant node by this node
					if (parent!=NULL){ // program statement list; do not replace
						h->n = th->n; //retarget pointers 
//						cout << "performed relinking" << endl;
					}
				}	
			}
		} else {
//			cout << "Entering recursion" << endl;
			makeComputationalTree( h->n, h->n->nodeArray, statementList);	
//			cout << "Leaving recursion" << endl;
		}
	h = h->next;
	}
}	

node* findVariableNodeByName(char* varName, nodeList* statementList){
	nodeList* sh = statementList;
	while (sh!=NULL){
		if ((strcmp(sh->n->name, varName)==0)) //if we found the node ...
			return sh->n; // return its pointer
		sh = sh->next;
	}
	//if we got here we didn't find the node, so there is a syntax error in the program
	cout<< "ERROR: output variable "<< varName << " was not declared;" << endl;
	exit(-1);	
}



bool isNodePresentInExpression(node* theNode, node* expression){
	nodeList* expressionOperandListHead = expression->nodeArray;
	cout << "is node " << theNode->name << " present in expression: "; 
	printExpression(expression); 
	cout << " ? " << endl;
	bool present = false;
	while (expressionOperandListHead!=NULL){
		if (expressionOperandListHead->n->nodeArray!=NULL)
			present = (present || isNodePresentInExpression(theNode, expressionOperandListHead->n));
		present = (present || (expressionOperandListHead->n == theNode));
		expressionOperandListHead = expressionOperandListHead->next;	
	}
	cout << "answer: " << present << endl;
	return present;
}



nodeList* createOuputList(nodeList* statementList, varList* variableList){
	/* for each output variable, find the corresponding node*/
	varList* vh = variableList;
	nodeList* useList;
	nodeList* potentialOutput = NULL;
	nodeList* properOutput = NULL;
		
	while (vh!=NULL){
		char* currentVariable = vh->name;
		cout << "processing variable " << vh->name << endl;	
		node* correspondingNode = findVariableNodeByName(currentVariable, statementList);	
		if (potentialOutput==NULL){
			correspondingNode->isOutput=true;
			potentialOutput = createNodeList(correspondingNode);
		}else{
			correspondingNode->isOutput=true;
			potentialOutput = appendNodeToListBack(correspondingNode, potentialOutput);
		}
		vh= vh->next;
	}
	cout << "Created Potential Output List" << endl;
	
	nodeList* poh = potentialOutput;
	nodeList* pohj = potentialOutput;
	/* iterate on the potential ouput nodes. For each one check if this node appears 
	in the calculation of some other node. If true, don't add it to the proper out list */
	while (poh!=NULL){
		bool isPresent=false;
		pohj = potentialOutput;
		while (pohj!=NULL){
			if (poh!=pohj)
				isPresent = (isPresent || isNodePresentInExpression(poh->n, pohj->n));
			pohj=pohj->next;
		}
		if (!isPresent){
			if (properOutput==NULL)	
				properOutput = createNodeList(poh->n);
			else
				properOutput = appendNodeToListBack( poh->n, properOutput);
		}
		poh = poh->next;
	}
	
	return properOutput;
}


