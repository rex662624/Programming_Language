%------------確認是不是質數(因為下面要找質數)

divisible(X,Y):-
N is Y*Y,N =< X,
X mod Y =:= 0.

divisible(X,Y):-
Y < X,Y1 is Y+1,
divisible(X,Y1).

is_prime(X):-
Y is 2, X > 1, 
\+divisible(X,Y).

%----------------------

% goldbach(N,S,L) : N是加起來的數字 S是小的數字 L是大的數字


goldbach(4,2,2).

%限定要偶數 而且>4
goldbach(N,S,L) :-   
N mod 2 =:= 0, 
N >4, 

%3表示從3這個質數找起，這樣下面nextprime就可以每次+2找下一個質數

goldbach(N,S,L,3).

%用剛剛傳進來的質數P(先假設它=S)測試，
%若N-S=L也是質數就表示找到兩個質數
%，若不是就要找下一個質數傳進來

goldbach(N,S,L,S) :- 
L is N - S, 
is_prime(L).

goldbach(N,S,L,P) :- 
P < N, next_prime(P,P1), 
goldbach(N,S,L,P1).

%P是目前的質數,P1是下一個質數，每次+2並確認找到的是不是質數

next_prime(P,P1) :- P1 is P + 2, is_prime(P1).
next_prime(P,P1) :- P2 is P + 2, next_prime(P2,P1).

%main
main:- readln(Input),member(X,Input),goldbach(X,Small,Large) , 
write(Small),write(" "),writeln(Large), halt.

:- initialization(main).
