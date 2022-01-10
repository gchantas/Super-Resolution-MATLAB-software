%Loads the array to GPU. Delete/comment line 5 to make this function the identity, so as to run on CPU if no GPU available.

function x=gpuArray2(x)

x=gpuArray(x)

