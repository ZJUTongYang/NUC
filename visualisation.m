mesh_tri = load("mesh_tri.txt");
mesh_ver = load("mesh_ver.txt");
result_path = load("result_path.txt");

mesh_tri = mesh_tri + 1;

plot3(result_path(:, 1), result_path(:, 2), result_path(:, 3), 'LineWidth', 3);