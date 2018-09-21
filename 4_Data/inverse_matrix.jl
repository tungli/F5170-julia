using LinearAlgebra

function myinv(mat::Matrix)
    s = size(mat)
    if s[1] != s[2]
        println("Not a square matrix!")
        return nothing
    elseif det(mat) == 0
        println("Matrix is singular!")
        return nothing
    else
        return inv(mat)
    end
end
