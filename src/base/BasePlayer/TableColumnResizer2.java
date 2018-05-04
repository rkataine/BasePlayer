package base.BasePlayer;

import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;

import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.event.MouseInputAdapter;
import javax.swing.table.TableColumn;



public class TableColumnResizer2 extends MouseInputAdapter
{ 
    public static Cursor resizeCursor = Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR); 
  
    private int mouseXOffset; 
    private Cursor otherCursor = resizeCursor; 
  
    private JTable table, table2; 
  
    public TableColumnResizer2(JTable table, JTable table2){ 
        this.table = table; 
        this.table2 = table2;
       
        table.addMouseListener(this); 
        table.addMouseMotionListener(this); 
    } 
  
    private boolean canResize(TableColumn column){ 
        return column != null
                && table.getTableHeader().getResizingAllowed() 
                && column.getResizable(); 
    } 
  
    private TableColumn getResizingColumn(JTable tablez, Point p){ 
        return getResizingColumn(table, p, table.columnAtPoint(p)); 
    } 
  
    private TableColumn getResizingColumn(JTable tablez, Point p, int column){ 
        if(column == -1){ 
            return null; 
        } 
        int row = table.rowAtPoint(p); 
        if(row==-1) 
            return null; 
        Rectangle r = table.getCellRect(row, column, true); 
        r.grow( -3, 0); 
        if(r.contains(p)) 
            return null; 
  
        int midPoint = r.x + r.width / 2; 
        int columnIndex; 
        if(table.getTableHeader().getComponentOrientation().isLeftToRight()) 
            columnIndex = (p.x < midPoint) ? column - 1 : column; 
        else
            columnIndex = (p.x < midPoint) ? column : column - 1; 
  
        if(columnIndex == -1) 
            return null;
        return tablez.getTableHeader().getColumnModel().getColumn(columnIndex);
    }
  
    public void mousePressed(MouseEvent e){ 
        table.getTableHeader().setDraggedColumn(null); 
        table.getTableHeader().setResizingColumn(null); 
        table.getTableHeader().setDraggedDistance(0); 
        
        Point p = e.getPoint(); 
  
        // First find which header cell was hit 
        int index = table.columnAtPoint(p); 
        if(index==-1) 
            return; 
  
        // The last 3 pixels + 3 pixels of next column are for resizing 
        TableColumn resizingColumn = getResizingColumn(table, p, index); 
        TableColumn resizingColumn2 = getResizingColumn(table2, p, index); 
        if(!canResize(resizingColumn)) 
            return; 
  
        table.getTableHeader().setResizingColumn(resizingColumn); 
        table2.getTableHeader().setResizingColumn(resizingColumn2);
        if(table.getTableHeader().getComponentOrientation().isLeftToRight()) 
            mouseXOffset = p.x - resizingColumn.getWidth(); 
        else
            mouseXOffset = p.x + resizingColumn.getWidth(); 
    } 
  
    private void swapCursor(){ 
        Cursor tmp = table.getCursor(); 
        table.setCursor(otherCursor); 
        otherCursor = tmp; 
    } 
  
    public void mouseMoved(MouseEvent e){ 
        if(canResize(getResizingColumn(table, e.getPoint())) 
           != (table.getCursor() == resizeCursor)){ 
            swapCursor(); 
        } 
    } 
  
    public void mouseDragged(MouseEvent e){ 
        int mouseX = e.getX(); 
  
        TableColumn resizingColumn = table.getTableHeader().getResizingColumn(); 
        
        TableColumn resizingColumn2 = table2.getTableHeader().getResizingColumn(); 
        boolean headerLeftToRight = table.getTableHeader().getComponentOrientation().isLeftToRight(); 
  
        if(resizingColumn != null){ 
            int oldWidth = resizingColumn.getWidth(); 
            int newWidth; 
            if(headerLeftToRight){ 
                newWidth = mouseX - mouseXOffset; 
            } else{
                newWidth = mouseXOffset - mouseX; 
            }
            resizingColumn.setWidth(newWidth); 
            resizingColumn2.setWidth(newWidth); 
            Container container; 
            if((table.getTableHeader().getParent() == null) || ((container = table.getTableHeader().getParent().getParent()) == null) || !(container instanceof JScrollPane)){ 
                return; 
            } 
  
            if(!container.getComponentOrientation().isLeftToRight() 
               && !headerLeftToRight){ 
                if(table != null){ 
                    JViewport viewport = ((JScrollPane)container).getViewport(); 
                    int viewportWidth = viewport.getWidth(); 
                    int diff = newWidth - oldWidth; 
                    int newHeaderWidth = table.getWidth() + diff; 
  
                    /* Resize a table */
                    Dimension tableSize = table.getSize(); 
                    tableSize.width += diff; 
                    table.setSize(tableSize); 
                    table2.setSize(tableSize); 
                    /* 
                     * If this table is in AUTO_RESIZE_OFF mode and has a horizontal 
                     * scrollbar, we need to update a view's position. 
                     */
                    if((newHeaderWidth >= viewportWidth) 
                       && (table.getAutoResizeMode() == JTable.AUTO_RESIZE_OFF)){ 
                        Point p = viewport.getViewPosition(); 
                        p.x = 
                                Math.max(0, Math.min(newHeaderWidth - viewportWidth, p.x + diff)); 
                        viewport.setViewPosition(p); 
  
                        /* Update the original X offset value. */
                        mouseXOffset += diff; 
                    } 
                } 
            } 
        } 
    } 
  
    public void mouseReleased(MouseEvent e){ 
        table.getTableHeader().setResizingColumn(null); 
        table.getTableHeader().setDraggedColumn(null); 
    } 
}