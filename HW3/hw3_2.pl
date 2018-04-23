
%relation_set 讀進n並設定relation
relationSet(N,Input):-%N是總數 input是每次讀進來的pair
	N>0,
	readln(Input),
	nth0(0,Input,P),	%取出input的第一個人當作是parents
	nth0(1,Input,C),
	assert(is_Parent(P,C)), %將這條規則加入fact
	N1 is N-1,
	relationSet(N1,_Input);N=0. %遇到N=0或是要繼續呼叫	

%getrelation
getRelation(T, P) :-
	T > 0,
	readln(P),
	nth0(0, P, A), nth0(1, P, B),
	lca(A,B),
	T0 is T - 1,
	getRelation(T0, _P);
	T = 0.

%ancestor

ancestor(A,C):-is_Parent(A,C).
ancestor(A,C):-is_Parent(P,C),ancestor(A,P).

lca(A,B) :- 
  A==B -> writeln(A);
  ancestor(A,B) -> writeln(A);
  is_Parent(X,A),lca(X,B).	

%mainfunction
main :-
	readln(N),
	member(Y,N),	
	relationSet(Y-1,_Input),%因為只有n-1個relation	
	
	%開始讀要求
	readln(T),
	member(Z, T),
	getRelation(Z, _A),
	halt.
:- initialization(main).
