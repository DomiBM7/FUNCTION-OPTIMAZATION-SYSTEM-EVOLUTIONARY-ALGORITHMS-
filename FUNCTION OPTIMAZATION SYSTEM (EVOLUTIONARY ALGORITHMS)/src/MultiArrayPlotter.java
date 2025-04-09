import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import javax.swing.JFrame;
import java.awt.Color;
import java.lang.reflect.Array;

public class MultiArrayPlotter extends JFrame {

    public MultiArrayPlotter(String title, double[]... arrays) {
        super(title);
        String[] labels = {"worstIndividuals", 
        		"averageIndividuals", 
        		"bestIndividuals", 
        		"rangeIndividuals",
        		"mutationIndividuals",
        			"selectivePressure",
        			"successRate",
        			"convergence",
        			"bestFitness",
        			"maxFitness",
        			"maxValues",
        		"deviation"};
        String[] labels2 = {};
        XYSeriesCollection dataset = new XYSeriesCollection();
        //represents a collection of dataset
        for (int i = 0; i < arrays.length; i++) {
        	XYSeries series = new XYSeries(labels[i]);
            for (int j = 0; j < arrays[i].length; j++) {
                series.add(j + 1, arrays[i][j]);
            }
            dataset.addSeries(series);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
            "Multi-Array Plot",
            "Generation",
            "Value",
            dataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );

        // Customize the chart and set colores
        chart.getXYPlot().setBackgroundPaint(Color.WHITE);
        chart.getXYPlot().setDomainGridlinePaint(Color.LIGHT_GRAY);
        chart.getXYPlot().setRangeGridlinePaint(Color.LIGHT_GRAY);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        setContentPane(chartPanel);
    }

    public static void plotArrays(String title, double[]... arrays) {
        MultiArrayPlotter plot = new MultiArrayPlotter(title, arrays);
        plot.pack();
        plot.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plot.setVisible(true);
    }

}