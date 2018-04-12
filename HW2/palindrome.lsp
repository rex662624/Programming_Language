(defun palindrome (l)
	
	(cond
			(	(equal l (reverse l))
				(progn
				(format t "~A is palindrome ~%" l)
				)		
			)
			(t
				(progn
				(format t "~A is not palindrome ~%" l) 
				)
			)
		)	
)

; MAIN FUNCTION
(let
	((list1(read)))
	;(setq list1  '(cat dog bird bird dog cat));
	(palindrome list1)
)
