using Base.Test
include("larStruct.jl")



@testset "box Tests" begin

	@testset "box Tests 2D" begin
	
		square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
		@test typeof(box(square))==Array{Array,1}
		@test length(box(square))==2
		@test length(box(square)[1])==2
	end
	
	
	@testset "box Tests 3D" begin
	
		cubes=([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]], [[[0], [1], [2], [3], [4], [5], [6], [7]], [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]], [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 4, 5], [2, 3, 6, 7], [0, 2, 4, 6], [1, 3, 5, 7]], [[0, 1, 2, 3, 4, 5, 6, 7]]])
		
		@test typeof(box(cubes))==Array{Array,1}
		@test length(box(cubes))==2
		@test length(box(cubes)[1])==3
	end
	
end




@testset "Struct Tests" begin

	square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
	@test Struct([square]).body==[square]
	@test Struct([square]).dim==length(square[1][1])
	@test Struct([square]).box==[[0,0],[1,1]]
	
end



@testset "embedStruct Tests" begin

	square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
	@test length(embedStruct(1)(Struct([square])).body[1][1][1])==length(Struct([square]).body[1][1][1])+1 
	#in questo caso n=1 in generale length(length(embedStruct(1)(x).body[1][1][1])=length(x.body[1][1][1])+n#
	@test length(embedStruct(3)(Struct([square])).body[1][1][1])==length(Struct([square]).body[1][1][1])+3
	@test typeof(embedStruct(1)(Struct([square])))==Struct	
	
end


@testset "embedTraversal Tests" begin

	square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
	@test length(embedTraversal(Struct([square]),Struct([square]),1,"New").body[1][1][1])==length(x.body[1][1][1])+1 
	#in questo caso n=1 in generale length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n#
	@test length(embedTraversal(Struct([square]),Struct([square]),3,"New").body[1][1][1])==length(x.body[1][1][1])+3
	@test typeof(embedTraversal(Struct([square]),Struct([square]),1,"New"))==Struct	
	
end


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
	 square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
	 table=larApply(t(-0.5,-0.5))(square)
	 structure=Struct([repeat([table,r(pi/2)],outer=2)...])
	 @test typeof(struct2lar(structure))==Tuple{Array{Any,1},Array{Array{Any,1},1}}
	 @test length(struct2lar(structure)[1][1])==2
	 end
	 
	 @testset "struct2lar 3D" begin
          BV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
          V=[[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
	 	  block=[V,BV]
	 	  structure=Struct(repeat([block,t(1,0,0)],outer=2));
	 	  @test typeof(struct2lar(structure))==Tuple{Array{Any,1},Array{Array{Any,1},1}}
	 	  @test length(struct2lar(structure)[1][1])==3

	 end

end



@testset "larRemoveVertices Tests" begin
	 V=[[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
	 FV=[[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 4, 5], [2, 3, 6, 7], [0, 2, 4, 6], [1, 3, 5, 7]]
	 @test typeof(larRemoveVertices(V,FV))==Tuple{Array{Any,1},Array{Any,1}}
	 @test length(larRemoveVertices(V,FV)[1])<= length(V)

end

	 
