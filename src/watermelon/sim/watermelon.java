package watermelon.sim;

// general utilities
import java.io.*;
import java.util.List;
import java.util.*;

import javax.tools.*;

import java.util.concurrent.*;
import java.net.URL;

// gui utility
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;

import javax.swing.*;

import watermelon.sim.seed;

public class watermelon {
	static String ROOT_DIR = "watermelon";

	// recompile .class file?
	static boolean recompile = true;

	// print more details?
	static boolean verbose = true;

	// Vary saturation of seeds by their score?
	static boolean enhancedColors = true;
	
	// Step by step trace
	static boolean trace = true;

	// enable gui
	static boolean gui = true;

	static double L;
	static double W;
	static ArrayList<Pair> treelist = new ArrayList<Pair>();
	static ArrayList<seed> seedlist = new ArrayList<seed>();

	static Player player;

	static double dimension = 100.0;
	static double distoseed = 1.999999999;
	static double distowall = 0.999999999;
	static double distotree = 1.999999999;
	static double s = 0.00;
	static double total = 0.00;

	static List<File> directoryFiles(String path, String extension) {
		List<File> allFiles = new ArrayList<File>();
		allFiles.add(new File(path));
		int index = 0;
		while (index != allFiles.size()) {
			File currentFile = allFiles.get(index);
			if (currentFile.isDirectory()) {
				allFiles.remove(index);
				for (File newFile : currentFile.listFiles())
					allFiles.add(newFile);
			} else if (!currentFile.getPath().endsWith(extension))
				allFiles.remove(index);
			else
				index++;
		}
		return allFiles;
	}

	// compile and load players dynamically
	//
	static Player loadPlayer(String group) {
		try {
			// get tools
			URL url = watermelon.class.getProtectionDomain().getCodeSource()
					.getLocation();
			// use the customized reloader, ensure clearing all static
			// information
			ClassLoader loader = new ClassReloader(url,
					watermelon.class.getClassLoader());
			if (loader == null)
				throw new Exception("Cannot load class loader");
			JavaCompiler compiler = null;
			StandardJavaFileManager fileManager = null;
			// get separator
			String sep = File.separator;
			// load players
			// search for compiled files
			File classFile = new File(ROOT_DIR + sep + group + sep
					+ "Player.class");
			System.err.println(classFile.getAbsolutePath());
			if (!classFile.exists() || recompile) {
				// delete all class files
				List<File> classFiles = directoryFiles(ROOT_DIR + sep + group,
						".class");
				System.err.print("Deleting " + classFiles.size()
						+ " class files...   ");
				for (File file : classFiles)
					file.delete();
				System.err.println("OK");
				if (compiler == null)
					compiler = ToolProvider.getSystemJavaCompiler();
				if (compiler == null)
					throw new Exception("Cannot load compiler");
				if (fileManager == null)
					fileManager = compiler.getStandardFileManager(null, null,
							null);
				if (fileManager == null)
					throw new Exception("Cannot load file manager");
				// compile all files
				List<File> javaFiles = directoryFiles(ROOT_DIR + sep + group,
						".java");
				System.err.print("Compiling " + javaFiles.size()
						+ " source files...   ");
				Iterable<? extends JavaFileObject> units = fileManager
						.getJavaFileObjectsFromFiles(javaFiles);
				boolean ok = compiler.getTask(null, fileManager, null, null,
						null, units).call();
				if (!ok)
					throw new Exception("Compile error");
				System.err.println("OK");
			}
			// load class
			System.err.print("Loading player class...   ");
			Class playerClass = loader.loadClass(ROOT_DIR + "." + group
					+ ".Player");
			System.err.println("OK");
			
			Player player = (Player) playerClass.newInstance();
			if (player == null)
				throw new Exception("Load error");
			else
				return player;

		} catch (Exception e) {
			e.printStackTrace(System.err);
			return null;
		}

	}

	// compute Euclidean distance between two points
	static double distance(seed a, Pair pair) {
		return Math.sqrt((a.x - pair.x) * (a.x - pair.x) + (a.y - pair.y) * (a.y - pair.y));
	}

	static double distanceseed(seed a, seed b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	/*
	 * static double vectorLength(double ox, double oy) {
	 * 
	 * return Math.sqrt(ox * ox + oy * oy); }
	 */
	void playgui() {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				OutpostUI ui = new OutpostUI();
				ui.createAndShowGUI();
			}
		});
	}

	class OutpostUI extends JPanel implements ActionListener {
		int FRAME_SIZE = 800;
		int FIELD_SIZE = 600;
		JFrame f;
		FieldPanel field;
		JButton next;
		//JButton next10;
		//JButton next50;
		JLabel label;
		JLabel label0;

		public OutpostUI() {
			setPreferredSize(new Dimension(FRAME_SIZE, FRAME_SIZE));
			setOpaque(false);
		}

		public void init() {
		}

		private boolean performOnce() {
			playStep();
			label.setText("Player Achieve the Score of " + total);
			label.setVisible(true);
			return true;
		}

		public void actionPerformed(ActionEvent e) {
			int steps = 0;

			if (e.getSource() == next)
				steps = 1;
			for (int i = 0; i < steps; ++i) {
				if (!performOnce())
					break;
			}

			repaint();
		}

		public void createAndShowGUI() {
			this.setLayout(null);

			f = new JFrame("Watermelon");
			field = new FieldPanel(1.0 * FIELD_SIZE / dimension);
			next = new JButton("Next");
			next.addActionListener(this);
			next.setBounds(0, 0, 100, 50);

			label = new JLabel();
			label.setVisible(false);
			label.setBounds(0, 60, 400, 50);
			label.setFont(new Font("Arial", Font.PLAIN, 15));

			field.setBounds(100, 100, FIELD_SIZE + 50, FIELD_SIZE + 50);


			this.add(next);
			this.add(label);
			this.add(field);

			f.add(this);

			f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			f.pack();
			f.setVisible(true);
		}

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
		}

	}

	class FieldPanel extends JPanel {
		double PSIZE = 10;
		double s;
		BasicStroke stroke = new BasicStroke(2.0f);
		double ox = 10.0;
		double oy = 10.0;

		public FieldPanel(double scale) {
			setOpaque(false);
			s = scale;
		}

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponent(g);

			Graphics2D g2 = (Graphics2D) g;
			g2.setStroke(stroke);
			double size = Math.max(W, L);
			// draw 2D rectangle
			g2.draw(new Rectangle2D.Double(ox, oy, dimension * s / size * W,
					dimension * s / size * L));
			
			double x_in = (dimension * s) / size;
			double y_in = (dimension * s) / size;
			
			for (int i=0; i<treelist.size(); i++) {
				g2.setPaint(Color.green);
				//g2.fill(new Rectangle2D.Double(ox+treelist.get(i).x,
					//	oy+treelist.get(i).y, x_in, y_in));
				drawTree(g2, treelist.get(i));
			}
			
			for (int i = 0; i < seedlist.size(); i++) {
				// drawPoint(g2, pointers[i]);
				drawPoint(g2, seedlist.get(i));
				// }
			}
		}
		
		public void drawTree(Graphics2D g2, Pair pr) {
				g2.setPaint(Color.GREEN);
			double size = Math.max(W, L);
			double x_in = (dimension * s) / size;
			double y_in = (dimension * s) / size;
			Ellipse2D e = new Ellipse2D.Double(ox + pr.x*x_in - x_in, oy + pr.y*y_in - y_in, 2*x_in , 2*y_in);
			g2.setStroke(stroke);
			g2.draw(e);
			g2.fill(e);

		}

		public void drawPoint(Graphics2D g2, seed sd) {
	        double saturation = sd.score*sd.score * 0.75 + 0.25;
	        
			if (sd.tetraploid == true) {
				if (enhancedColors)
					g2.setPaint(new Color(Color.HSBtoRGB((float) 0.6, (float) saturation, (float) 1.0)));
				else
					g2.setPaint(Color.BLACK);
			} else {
				if (enhancedColors)
					g2.setPaint(new Color(Color.HSBtoRGB((float) 0.0, (float) saturation, (float) 1.0)));
				else
					g2.setPaint(Color.MAGENTA);
			}
			
			double size = Math.max(W, L);
			double x_in = (dimension * s) / size;
			double y_in = (dimension * s) / size;
			
			Ellipse2D e = new Ellipse2D.Double(ox + sd.x*x_in - x_in, oy + sd.y*y_in - y_in, 2*x_in, 2*y_in);
			g2.setStroke(stroke);
			g2.draw(e);
			g2.fill(e);

		}

	}

	static boolean hascontained(ArrayList<Pair> list, Pair pr) {
		for (int j = 0; j < list.size(); j++) {

			if (list.get(j).equals(pr)) {
				return true;
			}
		}
		return false;
	}

	void playStep() {
		seedlist = player.move(treelist, W, L, s);
		if (validateseed(seedlist)) {
			System.out.printf("total score is %f\n", calculatescore());
		}

	}

	double calculatescore() {
		total = 0;
		
		for (int i = 0; i < seedlist.size(); i++) {
			double score;
			double chance = 0.0;
			double totaldis = 0.0;
			double difdis = 0.0;
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i) {
					totaldis = totaldis
							+ Math.pow(
									distanceseed(seedlist.get(i),
											seedlist.get(j)), -2);
				}
			}
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i
						&& ((seedlist.get(i).tetraploid && !seedlist.get(j).tetraploid) || (!seedlist
								.get(i).tetraploid && seedlist.get(j).tetraploid))) {
					difdis = difdis
							+ Math.pow(
									distanceseed(seedlist.get(i),
											seedlist.get(j)), -2);
				}
			}
			//System.out.println(totaldis);
			//System.out.println(difdis);
			chance = difdis / totaldis;
			score = chance + (1 - chance) * s;
			seedlist.get(i).score = score;
			total = total + score;
		}
		return total;
	}

	boolean validateseed(ArrayList<seed> seedlistin) {
		int nseeds = seedlistin.size();

		for (int i = 0; i < nseeds; i++) {
			for (int j = i + 1; j < nseeds; j++) {
				if (distanceseed(seedlistin.get(i), seedlistin.get(j)) < distoseed) {
					System.out.printf(
							"The distance between %d seed  %d seed is %f\n", i,
							j,
							distanceseed(seedlistin.get(i), seedlistin.get(j)));
					return false;
				}
			}
		}
		for (int i = 0; i < nseeds; i++) {
			if (seedlistin.get(i).x < 0 || seedlistin.get(i).x > W
					|| seedlistin.get(i).y < 0 || seedlistin.get(i).y > L) {
				System.out.printf(
						"The %d seed (%f, %f)  is out of the field\n", i,
						seedlistin.get(i).x, seedlistin.get(i).y);
				return false;
			}
			if (seedlistin.get(i).x < distowall
					|| W - seedlistin.get(i).x < distowall
					|| seedlistin.get(i).y < distowall
					|| L - seedlistin.get(i).y < distowall) {
				System.out.printf(
						"The %d seed (%f, %f) is too close to the wall\n", i,
						seedlistin.get(i).x, seedlistin.get(i).y);
				return false;
			}
		}
		for (int i = 0; i < treelist.size(); i++) {
			for (int j = 0; j < nseeds; j++) {
				if (distance(seedlistin.get(j), treelist.get(i)) < distotree) {
					System.out
							.printf("The %d seed (%f, %f) is too close to the tree (%f, %f), %f\n",
									j,
									seedlistin.get(j).x,
									seedlistin.get(j).y,
									treelist.get(i).x,
									treelist.get(i).y,
									distance(seedlistin.get(j), treelist.get(i)));
					return false;
				}
			}
		}
		return true;
	}

	void play() {
		/*
		 * while (tick <= MAX_TICKS) { if (endOfGame()) break; playStep(); }
		 * 
		 * if (tick > MAX_TICKS) { // Time out
		 * System.err.println("[ERROR] The player is time out!"); } else { //
		 * Achieve the goal
		 * System.err.println("[SUCCESS] The player achieves the goal in " +
		 * tick + " ticks."); }
		 */
	}

/*	void init() {
		grid = new Point[L * W];
		for (int i = 0; i < W; i++) {
			for (int j = 0; j < L; j++) {
				grid[i * L + j] = new Point(i, j, false);
			}
		}

	}*/

	watermelon() {
		this.player = player;
	}

	public void read(String map) {

		List<Pair> list = new ArrayList<Pair>();
		File file = new File(map);
		BufferedReader reader = null;
		int counter = 0;
		try {
			reader = new BufferedReader(new FileReader(file));
			String text = null;

			while ((text = reader.readLine()) != null) {
				java.util.List<String> item2 = new ArrayList();
				item2 = Arrays.asList(text.split(" "));
				ArrayList array_tmp0 = new ArrayList();
				Pair pr = new Pair();
				pr.x = Double.parseDouble(item2.get(0));
				pr.y = Double.parseDouble(item2.get(1));
				list.add(pr);
				counter = counter + 1;
				treelist.add(new Pair(pr));

			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException e) {
			}
		}
	}

	public static void main(String[] args) throws Exception {
		String map = null;
		String group = null;
		if (args.length > 0)
			map = args[0];
		if (args.length > 1)
			L = Double.parseDouble(args[1]);
		if (args.length > 2)
			W = Double.parseDouble(args[2]);
		if (args.length > 3)
			gui = Boolean.parseBoolean(args[3]);
		if (args.length > 4)
			group = args[4];
		if (args.length > 5)
			s = Double.parseDouble(args[5]);
		player = loadPlayer(group);
		watermelon game = new watermelon();
		//game.init();
		game.read(map);
		game.playgui();
	}
}