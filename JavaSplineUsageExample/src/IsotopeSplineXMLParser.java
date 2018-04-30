import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

public class IsotopeSplineXMLParser {

    Base64.Decoder decoder;

    public IsotopeSplineXMLParser() {
        decoder = Base64.getDecoder();
    }

    public Map<Integer, List<CubicSpline>> parse(String path) {

        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

        try {
            DocumentBuilder db = dbf.newDocumentBuilder();
            Document dom = db.parse(path);
            Map<Integer, List<CubicSpline>> models = parseDocument(dom);
            return models;
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    private Map<Integer, List<CubicSpline>> parseDocument(Document dom) {
        // get root element
        Element rootEle = dom.getDocumentElement();
        int numModels = Integer.parseInt(rootEle.getAttribute("maxIsotopeDepth"));
        Map<Integer, List<CubicSpline>> numSulfur2models = new HashMap();

        // get a nodelist of models
        NodeList nl = rootEle.getElementsByTagName("model");
        if (nl != null) {
            for (int i = 0; i < nl.getLength(); ++i) {

                int numSulfur = -1, isotope = 0;

                // get model element
                Element modelEle = (Element) nl.item(i);
                NamedNodeMap attributeMap = modelEle.getAttributes();

                for (int attIndex = 0; attIndex < attributeMap.getLength(); ++attIndex) {
                    if (attributeMap.item(attIndex).getNodeName() == "isotope")
                    {
                        isotope = Integer.parseInt(attributeMap.item(attIndex).getNodeValue());
                    } else if (attributeMap.item(attIndex).getNodeName() == "S") {
                        numSulfur = Integer.parseInt(attributeMap.item(attIndex).getNodeValue());
                    }
                }

                if (!numSulfur2models.containsKey(numSulfur)) {
                    numSulfur2models.put(numSulfur, Arrays.asList(new CubicSpline[numModels]));
                }

                numSulfur2models.get(numSulfur).set(isotope, parseModel(modelEle));
            }
        }

        return numSulfur2models;
    }

    private CubicSpline parseModel(Element modelEle) {
        double[] knots = decodeDoubleList(getTextValue(modelEle, "knots"));
        double[] coefficients = decodeDoubleList(getTextValue(modelEle, "coefficients"));

        List<Double> a = new ArrayList<Double>();
        List<Double> b = new ArrayList<Double>();
        List<Double> c = new ArrayList<Double>();
        List<Double> d = new ArrayList<Double>();
        List<Double> x = new ArrayList<Double>();

        for (int i = 0; i < coefficients.length; i+=4) {
            a.add(coefficients[i]);
            b.add(coefficients[i+1]);
            c.add(coefficients[i+2]);
            d.add(coefficients[i+3]);
        }

        for (int i = 0; i < knots.length; ++i) {
            x.add(knots[i]);
        }

        return new CubicSpline(a,b,c,d,x);
    }

    private double[] decodeDoubleList(String encoded) {

        ByteBuffer buff;
        byte[] decoded = decoder.decode(encoded.getBytes());

        double[] values = new double[decoded.length/8];

        for (int i = 0; i < decoded.length; i+=8) {
            byte[] data = new byte[8];
            for (int j = 0; j < 8; ++j) {
                data[j] = decoded[i+j];
            }

            buff = ByteBuffer.wrap(data);
            buff = buff.order(ByteOrder.LITTLE_ENDIAN);
            values[i/8] = buff.getDouble();
        }

        return values;
    }

    private String getTextValue(Element ele, String tagName) {
        String textVal = null;
        NodeList nl = ele.getElementsByTagName(tagName);
        if(nl != null && nl.getLength() > 0) {
            Element el = (Element)nl.item(0);
            textVal = el.getFirstChild().getNodeValue();
        }

        return textVal;
    }


}