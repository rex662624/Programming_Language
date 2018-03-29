(defun mergesort (numbers)
	 numbers
)




;main function
(let
;--------part1:declare variable-----
	 ((n (read))
	(numbers))
	
;--------part2:expression-----------
    ;set n element array
	(setf numbers(make-array n))

	;store number in array
	(do ((i 0 (+ i 1)))
		((>= i n))                   
		(setf (aref numbers i) (read)) 
    )

	;call merge sort
	
	;print sorted array
	(do ((i 0 (+ i 1)))
		((>= i n))
 		(format t "~A " (aref numbers i)) 
	)
	
	(format t "~%")
	
)
