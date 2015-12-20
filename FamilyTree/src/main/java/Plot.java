import com.lowagie.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;

import java.awt.*;
import java.io.FileNotFoundException;

/**
 * Created by hwestra on 10/20/15.
 */
public class Plot extends DefaultGraphics {


	private final FontMetrics metrics;
	DefaultTheme theme = new DefaultTheme();

	public Plot(String name, int width, int height) throws FileNotFoundException, DocumentException {
		super(name, width, height);
		g2d.setColor(theme.getDarkGrey());
		g2d.setFont(theme.getSmallFont());
		metrics = g2d.getFontMetrics();
	}

	public void drawInd(Individual ind, int y, int size) {
		Individual.Gender gender = ind.getGender();
		int x = ind.getXPixel();


		if (gender.equals(Individual.Gender.FEMALE)) {
			g2d.drawOval(x, y, size, size);
		} else {
			g2d.drawOval(x, y, size, size);
		}
		String lastname = ind.getLastName();
		String firstname = ind.getFirstName();
		int widthOfFirstName = metrics.stringWidth(firstname);
		int widthOfLastName = metrics.stringWidth(lastname);

		g2d.drawString(ind.getLastName(), x - (widthOfLastName / 2), y + size + 5);
		g2d.drawString(ind.getFirstName(), x - (widthOfFirstName / 2), y + size + 15);


	}
}
