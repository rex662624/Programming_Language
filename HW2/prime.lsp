; main functimn
(let 
	(
		(n (read))
		(i 0)
		(flag 0)
	)
	
	(do ((i 2 (+ i 1)))
		((> i (floor n 2)))
		(if (and (= flag 0) (= 0 (mod n i )))
			(progn
				(setf flag 1)
				
			)
			
		)
	)
	
		(cond
			(	(= flag 0 )
				(progn
				(format t "~A is a prime number ~%" n)
				)
					
			)
			
			(t
				(progn
				(format t "~A is not a prime number ~%" n) 

				)
			)

		)
	
)
