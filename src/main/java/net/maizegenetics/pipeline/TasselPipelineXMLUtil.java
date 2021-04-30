/*
 * TasselPipelineXMLUtil
 */
package net.maizegenetics.pipeline;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author Terry Casstevens
 */
public class TasselPipelineXMLUtil {

    private static final Logger myLogger = Logger.getLogger(TasselPipelineXMLUtil.class);

    private TasselPipelineXMLUtil() {
        // Utility Class
    }

    public static void writeArgsAsXML(String filename, String[] args) {

        try {
            DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder docBuilder = docFactory.newDocumentBuilder();

            Document doc = docBuilder.newDocument();
            Element rootElement = doc.createElement("TasselPipeline");
            doc.appendChild(rootElement);
            createXML(doc, rootElement, args);

            TransformerFactory transformerFactory = TransformerFactory.newInstance();
            Transformer transformer = transformerFactory.newTransformer();
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");
            transformer.setOutputProperty(OutputKeys.METHOD, "xml");
            DOMSource source = new DOMSource(doc);
            StreamResult result = new StreamResult(new File(filename));

            transformer.transform(source, result);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void createXML(Document doc, Element element, String[] args) throws IOException {

        int index = 0;

        while (index < args.length) {
            String current = args[index];
            if (!isFork(current)) {
                throw new IllegalArgumentException("TasselPipelineXMLUtil: createXML: this flag should be either -fork, -combine, or -runfork: " + current);
            }
            Element newElement = createTag(doc, element, current);

            while (true) {
                index++;
                if (index >= args.length) {
                    break;
                }
                if (isFork(args[index])) {
                    break;
                } else if (isSelfDescribingPlugin(args[index])) {
                    index = createSelfDescribingPluginXML(doc, newElement, args, index);
                } else if (isModifier(args[index])) {
                    index = createString(doc, newElement, args, index);
                } else {
                    index = createXML(doc, newElement, args, index);
                }
            }

        }

    }

    private static int createXML(Document doc, Element element, String[] args, int index) throws IOException {

        String current = args[index];
        Element newElement = createTag(doc, element, current);

        while (true) {
            index++;
            if (index >= args.length) {
                break;
            }
            if (isModifier(args[index])) {
                index = createString(doc, newElement, args, index);
            } else {
                break;
            }
        }

        return index - 1;

    }

    private static int createSelfDescribingPluginXML(Document doc, Element element, String[] args, int index) throws IOException {

        String current = args[index];
        Element newElement = createTag(doc, element, current);

        while (true) {
            index++;
            if (index >= args.length) {
                break;
            }
            if (args[index].equalsIgnoreCase("-endPlugin")) {
                index++;
                break;
            }
            if (args[index].startsWith("-runfork")) {
                break;
            }
            index = createString(doc, newElement, args, index);
        }

        return index - 1;

    }

    private static boolean isFork(String str) {

        if ((str.startsWith("-fork")) || (str.startsWith("-runfork")) || (str.startsWith("-combine"))) {
            return true;
        } else {
            return false;
        }

    }

    private static boolean isModifier(String str) {

        if (str.startsWith("-")) {

            TasselPipeline.FLAGS temp = null;
            try {
                temp = TasselPipeline.FLAGS.valueOf(str.substring(1));
            } catch (Exception e) {
                temp = null;
            }

            if ((str.startsWith("-fork")) || (str.startsWith("-runfork")) || (str.startsWith("-combine"))) {
                return false;
            } else if (temp != null) {
                return false;
            } else {
                return !isSelfDescribingPlugin(str);
            }

        } else {
            return true;
        }

    }

    private static boolean isSelfDescribingPlugin(String str) {

        if (str.startsWith("-")) {
            str = str.substring(1);
        }

        List<String> matches = Utils.getFullyQualifiedClassNames(str);
        for (String current : matches) {
            if (Plugin.isPlugin(current)) {
                return true;
            }
        }

        return false;

    }

    private static Element createTag(Document doc, Element element, String tag) {
        String str = tag.substring(tag.lastIndexOf('-') + 1);
        Element tagElement = doc.createElement(str);
        element.appendChild(tagElement);
        return tagElement;
    }

    private static int createString(Document doc, Element element, String[] args, int index) throws IOException {
        String current = args[index].substring(args[index].lastIndexOf('-') + 1);
        if (args[index].startsWith("-")) {
            Element newElement = createTag(doc, element, current);
            while (true) {
                index++;
                if (index >= args.length) {
                    return index - 1;
                }
                if (args[index].startsWith("-")) {
                    return index - 1;
                }
                newElement.appendChild(doc.createTextNode(args[index]));
            }
        } else {
            element.appendChild(doc.createTextNode(current));
            return index;
        }
    }

    public static String[][] readXMLAsArgs(String filename) {

        try {
            File fXmlFile = new File(filename);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            return readXMLAsArgs(doc);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TasselPipelineXMLUtil: readXMLAsArgs: Problem reading XML file: " + filename + "\n" + e.getMessage());
        }

    }

    public static String[][] readXMLAsArgsFromResource(String filename) {

        try {
            InputStream input = TasselPipelineXMLUtil.class.getResourceAsStream(filename);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(input);
            return readXMLAsArgs(doc);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TasselPipelineXMLUtil: readXMLAsArgsFromResource: Problem reading XML file: " + filename + "\n" + e.getMessage());
        }

    }

    private static String[][] readXMLAsArgs(Document doc) {

        try {

            doc.getDocumentElement().normalize();

            Element rootElement = doc.getDocumentElement();
            if (!(rootElement.getNodeName().equalsIgnoreCase("TasselPipeline"))) {
                throw new IllegalArgumentException("TasselPipelineXMLUtil: readXMLAsArgs: Root Node must be TasselPipeline: " + rootElement.getNodeName());
            }

            List<String> temp = new ArrayList<>();
            List<String> workflow = new ArrayList<>();
            NodeList children = rootElement.getChildNodes();
            String overallDescription = null;
            String citation = null;
            for (int i = 0; i < children.getLength(); i++) {
                Node current = children.item(i);
                if ((current.getNodeType() == Node.ELEMENT_NODE) && (current.getNodeName().trim().equals("workflow"))) {
                    NodeList descNodes = current.getChildNodes();
                    for (int j = 0; j < descNodes.getLength(); j++) {
                        Node currentDesc = descNodes.item(j);
                        if (currentDesc.getNodeType() == Node.TEXT_NODE) {
                            String value = currentDesc.getNodeValue().trim();
                            overallDescription = formatDescription(value);
                            break;
                        }
                    }
                } else if ((current.getNodeType() == Node.ELEMENT_NODE) && (current.getNodeName().trim().equals("citation"))) {
                    NodeList descNodes = current.getChildNodes();
                    for (int j = 0; j < descNodes.getLength(); j++) {
                        Node currentDesc = descNodes.item(j);
                        if (currentDesc.getNodeType() == Node.TEXT_NODE) {
                            String value = "Citation: " + currentDesc.getNodeValue().trim();
                            citation = formatDescription(value);
                            break;
                        }
                    }
                } else {
                    getFlags(current, temp, workflow);
                }
            }

            if (temp.size() != workflow.size()) {
                throw new IllegalStateException("TasselPipelineXMLUtil: readXMLAsArgs: number of flags and descriptions should be the same.");
            }

            String[] args = new String[temp.size()];
            temp.toArray(args);

            String[] descriptions = new String[workflow.size()];
            workflow.toArray(descriptions);

            String[][] result = new String[3][];
            result[0] = args;
            result[1] = descriptions;
            result[2] = new String[]{overallDescription, citation};
            return result;

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TasselPipelineXMLUtil: readXMLAsArgs: Problem converting XML to Args: " + e.getMessage());
        }

    }

    private static void getFlags(Node node, List<String> flags, List<String> workflow) {

        if (node.getNodeType() != Node.ELEMENT_NODE) {
            return;
        }

        String flagName = node.getNodeName().trim();
        final int workflowIndex = flags.size();
        flags.add("-" + flagName);
        workflow.add(null);
        NodeList children = node.getChildNodes();
        for (int i = 0; i < children.getLength(); i++) {
            Node current = children.item(i);
            if (current.getNodeType() == Node.TEXT_NODE) {
                String temp = current.getNodeValue().trim();
                if (temp.length() != 0) {
                    flags.add(temp);
                    workflow.add(null);
                }
            } else if ((current.getNodeType() == Node.ELEMENT_NODE) && (current.getNodeName().trim().equals("workflow"))) {
                NodeList descNodes = current.getChildNodes();
                for (int j = 0; j < descNodes.getLength(); j++) {
                    Node currentDesc = descNodes.item(j);
                    if (currentDesc.getNodeType() == Node.TEXT_NODE) {
                        workflow.remove(workflowIndex);
                        String value = currentDesc.getNodeValue().trim();
                        workflow.add(workflowIndex, formatDescription(value));
                        break;
                    }
                }
            } else {
                getFlags(current, flags, workflow);
            }
        }

        if (isSelfDescribingPlugin(flagName)) {
            flags.add("-endPlugin");
            workflow.add(null);
        }

    }

    private static final int DEFAULT_LINE_LENGTH = 50;

    private static String formatDescription(String description) {
        description = description.replaceAll("\\s+", " ");
        int count = 0;
        StringBuilder builder = new StringBuilder();
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > DEFAULT_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        return builder.toString();
    }
}
