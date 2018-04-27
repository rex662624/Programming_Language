
%relation_set 讀進n並設定relation
relationSet(N,Input):-%N是總數 input是每次讀進來的pair
	N>0,
	readln(Input),
	nth0(0,Input,P),	%取出input的第一個人當作是parents
	nth0(1,Input,C),
	assert(is_Parent(P,C)), %將P是c的parent這條規則加入fact
	N1 is N-1,
	relationSet(N1,_Input);N=0. %遇到N=0或是要繼續呼叫	

%getrelation 是一個for迴圈,根據使用者輸入的input去找答案
%因為最後一個輸入的答案會存在list最前面,所以最後要做reverse

getRelation(T,P,I1, I2) :-
	T > 0,
	readln(P),
	nth0(0, P, A), nth0(1, P, B),
	lca(A,B,I1,I2),
	T0 is T - 1,
	getRelation(T0, _P,I2,_S3);
	T = 0,
	reverse(I1,X,[]),%做reverse(因為最後輸入的會存在list的最前面)
	printlist(X).

%reverse 將list順序相反
reverse([],Z,Z).
reverse([H|T],Z,Acc):-reverse(T,Z,[H|Acc]).


%ancestor

ancestor(A,C):-is_Parent(A,C).
ancestor(A,C):-is_Parent(P,C),ancestor(A,P).

%lowest_commom_ancestor
lca(A,B,I1,I2) :- 
  A==B -> append([A],I1,I2);
  ancestor(A,B) -> append([A],I1,I2);
  is_Parent(X,A),lca(X,B,I1,I2).	

%將List print出來

printlist([]).
printlist([X|List]) :-
    writeln(X),
    printlist(List).

%mainfunction
main :-
	readln(N),
	member(Y,N),	
	relationSet(Y-1,_Input),%因為只有n-1個relation	
	
	%開始讀要求
	readln(T),
	member(Z, T),
	getRelation(Z, _A,_I1,_I2),
	halt.

:- initialization(main).
