using Base.Test
@testset "removeDups Tests" begin

	 @testset "removeDups 3D" begin
	 	   CW=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],[4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
		    @test length(removeDups(CW))<= length(CW)
		     @test typeof(removeDups(CW))==Array{Array{Int64,1},1}
	  end
	  
	   @testset "removeDups 2D" begin
	   	    CW=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]
		    @test length(removeDups(CW))<= length(CW)
		     @test typeof(removeDups(CW))==Array{Array{Int64,1},1}
	    end
end


@testset "struct2lar" begin

	 @testset "struct2lar 2D" begin


	 end
	 
	 @testset "struct2lar 3D" begin
	 	  block
	 	  structure=S

	 end

end





	 