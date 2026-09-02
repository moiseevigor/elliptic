function x = gather(g)
  if isa(g, 'gpuArray'), x = g.d; else, error('gather: argument is not a device array (ocl errors here too)'); end
end
