
% 宣告事實（哪些node相連）
%因為是雙向關係所以輸入一次要加入兩個fact ba 與ab
relationSet(N, Edge) :-
	N > 0,M is N-1,
	readln(Edge),
	nth0(0, Edge, A), nth0(1, Edge, B),
	assert(has(A, B)),
	assert(has(B, A)),
	relationSet(M, _P);N = 0.

% getRelation. for迴圈輸入值,與作業3.2動作相同
getRelation(T, L, P,I1,I2) :-
	T > 0,T0 is T - 1,
	readln(P),
	nth0(0, P, A), nth0(1, P, B),
	path(A,B,I1,I2),
	getRelation(T0, L, _A,I2,_S3);
	T = 0,
	reverse(I1,X,[]),%做reverse(因為最後輸入的會存在list的最前面)
	printlist(X).

%reverse 
reverse([],Z,Z).
reverse([H|T],Z,Acc):-reverse(T,Z,[H|Acc]).

%reachable
%建立雙向關係
%bi-has(X,Y):-has(X,Y).
%bi-has(Y,X):-has(Y,X).

path(A,B,I1,I2):-
	%要嘛走得到 要嘛走不到
	walk(A,B,[],I1,I2);append(['No'],I1,I2).%writeln('NO').

walk(A,B,V,I1,I2):-
	has(A,X),
	not(member(X,V)),
	(
	B=X -> append(['Yes'],I1,I2);%->writeln('Yes');
	B=\=X ->walk(X,B,[A|V],I1,I2)
	).
	

%將List print出來

printlist([]).
printlist([X|List]) :-
    writeln(X),
    printlist(List).


%mainfunction
main :-
	readln(L),
	nth0(0, L, N), nth0(1, L, E),	
	relationSet(E,_List),
	
	%開始讀要求
	readln(T),
	member(Z, T),
	getRelation(Z,N,_I1,_I2,_S3),
	halt.

:- initialization(main).
